
#include "../include/readtiff.h"
#include "mfem.hpp"

Constraints::Constraints() : Row_begin(0), Row_end(-1), Column_begin(0), Column_end(-1), Depth_begin(0), Depth_end(-1) {}

Constraints::Constraints(int row0, int row1, int col0, int col1, int depth0, int depth1)
    : Row_begin(row0), Row_end(row1), Column_begin(col0), Column_end(col1), Depth_begin(depth0), Depth_end(depth1) {}

TIFFReader::TIFFReader(const char* filePath, const Constraints& constraints) {
    this->filePath = filePath;
    tiff = TIFFOpen(filePath, "r");
    if (tiff == nullptr) {
        std::cerr << "Could not open TIFF file: " << filePath << std::endl;
        exit(1);
    }
    calculateNumPages();
    readTIFFFields();
    setConstraints(constraints);
    TIFFSetDirectory(tiff, 0);

    uint16 spp=0, bps=0, photo=0, planar=0;
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photo);
    TIFFGetField(tiff, TIFFTAG_PLANARCONFIG, &planar);

    std::cout << "spp=" << spp
            << " bps=" << bps
            << " photometric=" << photo
            << " planar=" << planar << "\n";

}



void TIFFReader::readinfo() {

    imageData.resize(constraints.Depth_end - constraints.Depth_begin);
    for (int page = constraints.Depth_begin; page < constraints.Depth_end; ++page) {
        imageData[page - constraints.Depth_begin].resize(constraints.Row_end - constraints.Row_begin);
        for (int row = constraints.Row_begin; row < constraints.Row_end; ++row) {
            imageData[page - constraints.Depth_begin][row - constraints.Row_begin].resize(constraints.Column_end - constraints.Column_begin);
        }
    }

    if (mfem::Mpi::WorldRank() == 0) { std::cout << "Original TIFF Info - Width: " << Width << ", Height: " << Height
              << ", NumPages: " << numPages << std::endl;}

    for (int page = 0; page < numPages; page++) {
        if (!(page > constraints.Depth_begin - 1 && page < constraints.Depth_end)) {
            TIFFSetDirectory(tiff, page);
            continue;
        }

        // if (mfem::Mpi::WorldRank() == 0) {std::cerr << "page: " << page << ", READING\n";}
        TIFFSetDirectory(tiff, page);

        // --- Read per-page photometric + spp (single slice is RGBA, stack is grayscale) ---
        uint16 photo = 0;
        TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photo);

        uint16 spp = 1;
        TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &spp);

        tdata_t buf = _TIFFmalloc(TIFFScanlineSize(tiff));

        // IMPORTANT: read scanline once per ROW (not once per (row,col))
        for (int row = constraints.Row_begin; row < constraints.Row_end; row++) {
            TIFFReadScanline(tiff, buf, row);
            uint8* p = (uint8*)buf; // scanline bytes

            for (int col = constraints.Column_begin; col < constraints.Column_end; col++) {

                // --- Correctly decode pixel value for grayscale vs RGBA ---
                // uint8 gray;
                // if (spp == 1) {
                //     gray = p[col];
                // } else {
                //     // packed RGBA: [R G B A] [R G B A] ...
                //     const int idx = (int)spp * col;
                //     const uint8 r = p[idx + 0];
                //     const uint8 g = p[idx + 1];
                //     const uint8 b = p[idx + 2];
                //     gray = static_cast<uint8>(0.299*r + 0.587*g + 0.114*b);
                // }

                uint8 gray;

                if (spp == 1) {
                    // grayscale
                    gray = p[col];
                }
                else if (spp == 2) {
                    // grayscale + alpha: [G A] [G A] ...
                    gray = p[2*col];      // ignore alpha
                }
                else {
                    // RGB/RGBA: [R G B (A)] ...
                    const int idx = (int)spp * col;
                    const uint8 r = p[idx + 0];
                    const uint8 g = p[idx + 1];
                    const uint8 b = p[idx + 2];
                    gray = static_cast<uint8>(0.299*r + 0.587*g + 0.114*b);
                }

                int solid;
                if (photo == PHOTOMETRIC_MINISWHITE) {
                    // 0=white, 255=black  -> solid is WHITE, so solid when gray is LOW
                    solid = (gray <= 127) ? 1 : 0;
                } else { // PHOTOMETRIC_MINISBLACK (and most other cases)
                    // 0=black, 255=white  -> solid is WHITE, so solid when gray is HIGH
                    solid = (gray >= 127) ? 1 : 0;
                }

                imageData[page - constraints.Depth_begin]
                         [row  - constraints.Row_begin]
                         [col  - constraints.Column_begin] = solid;
            }
        }

        _TIFFfree(buf);

        if (mfem::Mpi::WorldRank() == 0) {std::cout << "[TIFFReader] Constrained dimensions:\n"
          << "  Pages   : " << (constraints.Depth_end  - constraints.Depth_begin) << "\n"
          << "  Rows    : " << (constraints.Row_end    - constraints.Row_begin)   << "\n"
          << "  Columns : " << (constraints.Column_end - constraints.Column_begin) << std::endl;}

    }
}

TIFFReader::~TIFFReader() {
    if (tiff != nullptr) {
        TIFFClose(tiff);
    }
}

void TIFFReader::calculateNumPages() {
    numPages = 0;
    do {
        numPages++;
    } while (TIFFReadDirectory(tiff));
}

void TIFFReader::readTIFFFields() {
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &Width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &Height);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
}

void TIFFReader::setConstraints(const Constraints& constraints) {        
    this->constraints.Row_begin = std::max(constraints.Row_begin, 0);
    this->constraints.Column_begin = std::max(constraints.Column_begin, 0);
    this->constraints.Depth_begin = std::max(constraints.Depth_begin, 0);

    // NOTE: keeping your structure, but fixing the Width/Height swap:
    // Rows should clamp to Height, Columns should clamp to Width.

    if (constraints.Row_end > Height) {
        std::cerr << "Row_end (" << constraints.Row_end << ") exceeds Height (" << Height << "). Setting to Height.\n";
        this->constraints.Row_end = Height;
    } else if (constraints.Row_end < this->constraints.Row_begin) {
        std::cerr << "Row_end (" << constraints.Row_end << ") is less than Row_begin (" << this->constraints.Row_begin << "). Setting to Height.\n";
        this->constraints.Row_end = Height;
    } else {
        this->constraints.Row_end = (constraints.Row_end > -1) ? constraints.Row_end : Height;
    }

    if (constraints.Column_end > Width) {
        std::cerr << "Column_end (" << constraints.Column_end << ") exceeds Width (" << Width << "). Setting to Width.\n";
        this->constraints.Column_end = Width;
    } else if (constraints.Column_end < this->constraints.Column_begin) {
        std::cerr << "Column_end (" << constraints.Column_end << ") is less than Column_begin (" << this->constraints.Column_begin << "). Setting to Width.\n";
        this->constraints.Column_end = Width;
    } else {
        this->constraints.Column_end = (constraints.Column_end > -1) ? constraints.Column_end : Width;
    }

    if (constraints.Depth_end > numPages) {
        std::cerr << "Depth_end (" << constraints.Depth_end << ") exceeds total depth (" << numPages << "). Setting to number of pages.\n";
        this->constraints.Depth_end = numPages;
    } else if (constraints.Depth_end < this->constraints.Depth_begin) {
        std::cerr << "Depth_end (" << constraints.Depth_end << ") is less than Depth_begin (" << this->constraints.Depth_begin << "). Setting to total depth.\n";
        this->constraints.Depth_end = numPages;
    } else {
        this->constraints.Depth_end = (constraints.Depth_end > -1) ? constraints.Depth_end : numPages;
    }
    
}

std::vector<std::vector<std::vector<int>>> TIFFReader::getImageData() {
    return imageData;
}
