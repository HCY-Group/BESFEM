#include "readtiff.h"

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
}

void TIFFReader::readinfo() {

    imageData.resize(constraints.Depth_end - constraints.Depth_begin);
    for (int page = constraints.Depth_begin; page < constraints.Depth_end; ++page) {
        imageData[page - constraints.Depth_begin].resize(constraints.Row_end - constraints.Row_begin);
        for (int row = constraints.Row_begin; row < constraints.Row_end; ++row) {
            imageData[page - constraints.Depth_begin][row - constraints.Row_begin].resize(constraints.Column_end - constraints.Column_begin);
        }
    }

    for (int page = 0; page < numPages - 1; page++) {
        if (!(page > constraints.Depth_begin - 1 && page < constraints.Depth_end)) {
            std::cerr << "page: " << page + 1 << ", SKIPPING\n";
            TIFFSetDirectory(tiff, page);
            continue;
        }

        std::cerr << "page: " << page << ", READING\n";
        TIFFSetDirectory(tiff, page);
        tdata_t buf = _TIFFmalloc(TIFFScanlineSize(tiff));
        for (int col = constraints.Column_begin; col < constraints.Column_end; col++) {
            for (int row = constraints.Row_begin; row < constraints.Row_end; row++) {
                TIFFReadScanline(tiff, buf, row);
                uint8* pixel = (uint8*)buf; //cast void pointer
                imageData[page - constraints.Depth_begin][row - constraints.Row_begin][col - constraints.Column_begin] = 1 - static_cast<int>(pixel[col]) / 255;
            }
        }
        _TIFFfree(buf);
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

    if (constraints.Row_end > Width) {
        std::cerr << "Row_end (" << constraints.Row_end << ") exceeds Width (" << Width << "). Setting to Width.\n";
        this->constraints.Row_end = Width;
    } else if (constraints.Row_end < this->constraints.Row_begin) {
        std::cerr << "Row_end (" << constraints.Row_end << ") is less than Row_begin (" << this->constraints.Row_begin << "). Setting to Width.\n";
        this->constraints.Row_end = Width;
    } else {
        this->constraints.Row_end = (constraints.Row_end > -1) ? constraints.Row_end : Width;
    }

    if (constraints.Column_end > Height) {
        std::cerr << "Column_end (" << constraints.Column_end << ") exceeds Height (" << Height << "). Setting to Height.\n";
        this->constraints.Column_end = Height;
    } else if (constraints.Column_end < this->constraints.Column_begin) {
        std::cerr << "Column_end (" << constraints.Column_end << ") is less than Column_begin (" << this->constraints.Column_begin << "). Setting to Height.\n";
        this->constraints.Column_end = Height;
    } else {
        this->constraints.Column_end = (constraints.Column_end > -1) ? constraints.Column_end : Height;
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
       
/*
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <tiff_file> [<Row_begin> <Row_end> <Column_begin> <Column_end> <Depth_begin> <Depth_end>]\n";
        return 1;
    }
    Constraints args;
    for (int i = 2; i < argc; ++i) {
        int value = std::stoi(argv[i]);        switch (i) {
            case 2:
                args.Row_begin = value;
                break;
            case 3:
                args.Row_end = value;
                break;
            case 4:
                args.Column_begin = value;
                break;
            case 5:
                args.Column_end = value;
                break;
            case 6:
                args.Depth_begin = value;
                break;
            case 7:
                args.Depth_end = value;
                break;
            default:
                break;
        }
    }

    TIFFReader reader(argv[1], args);
    reader.readinfo();
    reader.getImageData();
    return 0;
}
*/
