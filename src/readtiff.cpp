
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

    uint16 spp = 0;
    uint16 bps = 0;
    uint16 photo = 0;
    uint16 planar = 0;
    uint16 sample_format = SAMPLEFORMAT_UINT;

    TIFFGetFieldDefaulted(
        tiff,
        TIFFTAG_SAMPLEFORMAT,
        &sample_format);

    if (mfem::Mpi::WorldRank() == 0)
    {
        std::cout
            << "spp=" << spp
            << " bps=" << bps
            << " sampleformat=" << sample_format
            << " photometric=" << photo
            << " planar=" << planar
            << "\n";
    }
    

    // TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &spp);
    // TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bps);
    // TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photo);
    // TIFFGetField(tiff, TIFFTAG_PLANARCONFIG, &planar);

    // std::cout << "spp=" << spp
    //         << " bps=" << bps
    //         << " photometric=" << photo
    //         << " planar=" << planar << "\n";

}

void TIFFReader::readinfo()
{
    const int nz = constraints.Depth_end  - constraints.Depth_begin;
    const int ny = constraints.Row_end    - constraints.Row_begin;
    const int nx = constraints.Column_end - constraints.Column_begin;

    imageData.assign(
        nz,
        std::vector<std::vector<int>>(
            ny,
            std::vector<int>(nx, 0)
        )
    );

    if (mfem::Mpi::WorldRank() == 0) {
        std::cout << "Original TIFF Info - Width: " << Width
                  << ", Height: " << Height
                  << ", NumPages: " << numPages << std::endl;
    }

    std::set<int> observed_values;

    uint16 first_photo = 0;
    uint16 first_spp   = 1;
    bool metadata_set  = false;

    for (int page = 0; page < numPages; page++) {
        if (!(page > constraints.Depth_begin - 1 && page < constraints.Depth_end)) {
            TIFFSetDirectory(tiff, page);
            continue;
        }

        TIFFSetDirectory(tiff, page);

        uint16 photo = 0;
        TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photo);

        uint16 spp = 1;
        TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &spp);

        if (!metadata_set) {
            first_photo = photo;
            first_spp   = spp;
            metadata_set = true;
        }

        uint16 bits_per_sample = 8;
        uint16 sample_format = SAMPLEFORMAT_UINT;

        TIFFGetFieldDefaulted(
            tiff,
            TIFFTAG_BITSPERSAMPLE,
            &bits_per_sample);

        TIFFGetFieldDefaulted(
            tiff,
            TIFFTAG_SAMPLEFORMAT,
            &sample_format);

        tdata_t buf =
            _TIFFmalloc(
                TIFFScanlineSize(tiff));

        if (buf == nullptr)
        {
            throw std::runtime_error(
                "Could not allocate TIFF scanline buffer.");
        }

        for (int row = constraints.Row_begin;
            row < constraints.Row_end;
            ++row)
        {
            if (TIFFReadScanline(
                    tiff,
                    buf,
                    row,
                    0) < 0)
            {
                _TIFFfree(buf);

                throw std::runtime_error(
                    "Failed to read TIFF scanline.");
            }

            for (int col = constraints.Column_begin;
                col < constraints.Column_end;
                ++col)
            {
                int value = 0;

                // -------------------------------------------------
                // Single-channel 8-bit TIFF
                // -------------------------------------------------
                if (spp == 1 &&
                    bits_per_sample == 8)
                {
                    const uint8_t *pixels =
                        static_cast<const uint8_t *>(buf);

                    value =
                        static_cast<int>(
                            pixels[col]);
                }

                // -------------------------------------------------
                // Single-channel signed 16-bit TIFF
                // -------------------------------------------------
                else if (
                    spp == 1 &&
                    bits_per_sample == 16 &&
                    sample_format == SAMPLEFORMAT_INT)
                {
                    const int16_t *pixels =
                        static_cast<const int16_t *>(buf);

                    value =
                        static_cast<int>(
                            pixels[col]);
                }

                // -------------------------------------------------
                // Single-channel unsigned 16-bit TIFF
                // -------------------------------------------------
                else if (
                    spp == 1 &&
                    bits_per_sample == 16 &&
                    sample_format == SAMPLEFORMAT_UINT)
                {
                    const uint16_t *pixels =
                        static_cast<const uint16_t *>(buf);

                    value =
                        static_cast<int>(
                            pixels[col]);
                }

                // -------------------------------------------------
                // RGB 8-bit TIFF
                // -------------------------------------------------
                else if (
                    spp >= 3 &&
                    bits_per_sample == 8)
                {
                    const uint8_t *pixels =
                        static_cast<const uint8_t *>(buf);

                    const int idx =
                        static_cast<int>(spp) * col;

                    const uint8_t r =
                        pixels[idx + 0];

                    const uint8_t g =
                        pixels[idx + 1];

                    const uint8_t b =
                        pixels[idx + 2];

                    const bool is_white =
                        r > 240 &&
                        g > 240 &&
                        b > 240;

                    value =
                        is_white ? 0 : 1;
                }
                else
                {
                    _TIFFfree(buf);

                    std::stringstream error;

                    error
                        << "Unsupported TIFF format: "
                        << "samples per pixel = "
                        << spp
                        << ", bits per sample = "
                        << bits_per_sample
                        << ", sample format = "
                        << sample_format
                        << ".";

                    throw std::runtime_error(
                        error.str());
                }

                imageData
                    [page - constraints.Depth_begin]
                    [row - constraints.Row_begin]
                    [col - constraints.Column_begin] =
                        value;

                observed_values.insert(
                    value);
            }
        }

        _TIFFfree(buf);
    }

    const int min_value = *observed_values.begin();
    const int max_value = *observed_values.rbegin();

    const bool is_binary_01 =
        observed_values.size() <= 2 &&
        min_value == 0 &&
        max_value == 1;

    const bool is_binary_255 =
        observed_values.size() <= 2 &&
        min_value == 0 &&
        max_value == 255;

    const bool is_grayscale =
        observed_values.size() > 20 &&
        min_value == 0 &&
        max_value <= 255;

    const bool has_negative_labels =
        min_value < 0;

    const bool has_positive_labels =
        max_value > 1;

    const bool is_label_tiff =
        has_negative_labels ||
        (
            !is_binary_01 &&
            !is_binary_255 &&
            !is_grayscale &&
            has_positive_labels
        );

    if (is_label_tiff) {
        if (mfem::Mpi::WorldRank() == 0) {
            std::cout << "[TIFFReader] Detected label TIFF. Keeping labels directly.\n";
        }
    }
    else if (is_binary_01) {
        if (mfem::Mpi::WorldRank() == 0) {
            std::cout << "[TIFFReader] Detected binary 0/1 TIFF. Keeping as 0/1.\n";
        }
    }
    else {
        if (mfem::Mpi::WorldRank() == 0) {
            std::cout << "[TIFFReader] Detected binary/grayscale TIFF. Converting to solid mask.\n";
        }

        for (auto &slice : imageData) {
            for (auto &row : slice) {
                for (int &v : row) {

                    if (first_photo == PHOTOMETRIC_MINISBLACK) {
                        // MINISBLACK: black is low value
                        // so black particle means v < 127
                        v = (v < 127) ? 1 : 0;
                    }
                    else if (first_photo == PHOTOMETRIC_MINISWHITE) {
                        // MINISWHITE: black is high value
                        // so black particle means v > 127
                        v = (v < 127) ? 1 : 0;
                    }
                    else {
                        v = (v < 127) ? 1 : 0;
                    }
                }
            }
        }
    }

    if (mfem::Mpi::WorldRank() == 0) {
        std::cout << "[TIFFReader] Constrained dimensions:\n"
                  << "  Pages   : " << nz << "\n"
                  << "  Rows    : " << ny << "\n"
                  << "  Columns : " << nx << "\n"
                  << "  Values found: ";

        for (int v : observed_values) std::cout << v << " ";
        std::cout << std::endl;
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
