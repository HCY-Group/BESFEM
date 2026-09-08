#include <tiffio.h>
#include <iostream>
#include <vector>
#include "SimulationConfig.hpp"


struct Constraints {
    int Row_begin, Row_end, Column_begin, Column_end, Depth_begin, Depth_end;
    Constraints();
    Constraints(int row0, int row1, int col0, int col1, int depth0, int depth1);
};

class TIFFReader {
public:
    const SimulationConfig& cfg;

    TIFFReader(const char* filePath, const Constraints& constraints, const SimulationConfig& cfg);
    void readinfo();
    ~TIFFReader();
    void calculateNumPages();
    void readTIFFFields();
    void setConstraints(const Constraints& constraints);
    std::vector<std::vector<std::vector<int>>> getImageData();

private:
    std::vector<std::vector<std::vector<int>>> imageData;
    const char* filePath;
    TIFF* tiff;
    int numPages;
    int Width, Height, samplesPerPixel, bitsPerSample;
    Constraints constraints;

};