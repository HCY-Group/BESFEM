/**
 * @file readtiff.h
 * @brief Read cropped TIFF stacks as particle labels or solid-phase masks.
 */
#include <tiffio.h>
#include <iostream>
#include <vector>
#include "SimulationConfig.hpp"

/**
 * @brief Zero-based crop bounds for rows, columns, and TIFF pages.
 *
 * Begin indices are inclusive and end indices are exclusive. An end value
 * of -1 requests the corresponding full dimension when applied by TIFFReader.
 */
struct Constraints {
    int Row_begin, ///< Inclusive first row (y).
        Row_end, ///< Exclusive last row; -1 requests the image height.
        Column_begin, ///< Inclusive first column (x).
        Column_end, ///< Exclusive last column; -1 requests the image width.
        Depth_begin, ///< Inclusive first TIFF page (z).
        Depth_end; ///< Exclusive last page; -1 requests the page count.

    /** @brief Set all begin indices to zero and all end indices to -1. */
    Constraints();

    /**
     * @brief Store explicit crop bounds without validating them.
     * @param row0 Inclusive first row.
     * @param row1 Exclusive last row, or -1 for the full height.
     * @param col0 Inclusive first column.
     * @param col1 Exclusive last column, or -1 for the full width.
     * @param depth0 Inclusive first page.
     * @param depth1 Exclusive last page, or -1 for the full stack depth.
     */
    Constraints(int row0, int row1, int col0, int col1, int depth0, int depth1);
};

/**
 * @brief Owns a TIFF handle and reads a cropped stack into integer voxel data.
 *
 * Construction opens the file and loads metadata. Call readinfo() to populate
 * the voxel data, then getImageData() to obtain a copy indexed by [z][y][x].
 * Detected label images retain their values; other images may be converted
 * to a binary solid mask using the configuration and TIFF photometric tags.
 */
class TIFFReader {
public:
    const SimulationConfig& cfg; ///< Borrowed configuration; must outlive this reader.

    /**
     * @brief Open the TIFF, count pages, read metadata, and apply crop bounds.
     * @param filePath Path to the TIFF file; the pointer is stored without copying.
     * @param constraints Requested crop bounds.
     * @param cfg Configuration specifying particle-color interpretation.
     * @note Exits the process with status 1 if TIFFOpen fails.
     * @note Pixel data is not loaded until readinfo() is called.
     */
    TIFFReader(const char* filePath, const Constraints& constraints, const SimulationConfig& cfg);

    /**
     * @brief Read the selected region and classify or convert its pixel values.
     *
     * Reads single-channel 8-bit values as unsigned, single-channel signed or
     * unsigned 16-bit values, and interleaved 8-bit RGB data with at least three
     * samples per pixel. RGB pixels with all three channels above 240 become
     * zero; other RGB pixels become one before image classification.
     *
     * A value-based heuristic identifies label images and binary 0/1 images
     * to preserve. Remaining images are thresholded at 127 using photometric
     * interpretation and cfg.particle_color. Replaces previously loaded data.
     * @pre Crop bounds select a nonempty region within the image dimensions.
     * @pre Selected pages have compatible dimensions and scanline layouts.
     * @throws std::runtime_error If scanline allocation or reading fails, or
     *         a sample format is unsupported by the decoding branches.
     */
    void readinfo();

    /** @brief Close the owned TIFF handle when it is non-null. */
    ~TIFFReader();

    /**
     * @brief Count directories from the current TIFF directory to the last.
     * @note Advances the TIFF directory cursor and updates numPages. To count
     *       the entire stack, the current directory must be the first page.
     */
    void calculateNumPages();

    /**
     * @brief Read width, height, samples per pixel, and bit depth from the current page.
     * @note Updates cached metadata without loading pixel data.
     */
    void readTIFFFields();

    /**
     * @brief Normalize and store crop bounds using cached image dimensions.
     *
     * Negative begins become zero. Ends beyond the dimension or below the
     * normalized begin become the full dimension, including the -1 sentinel.
     * Begin indices are not capped at the image dimensions, and empty ranges
     * are not rejected. Does not reload imageData.
     * @param constraints Requested row, column, and page bounds.
     */
    void setConstraints(const Constraints& constraints);

    /**
     * @brief Return a copy of the last loaded cropped voxel data.
     * @return Integer values indexed by [page][row][column], relative to the
     *         crop origin; empty if readinfo() has not yet populated the data.
     */
    std::vector<std::vector<std::vector<int>>> getImageData();

private:
    std::vector<std::vector<std::vector<int>>> imageData; ///< Cropped labels or mask, indexed by [z][y][x].
    const char* filePath; ///< Borrowed file-path pointer supplied to the constructor.
    TIFF* tiff; ///< Owned libtiff handle, closed by the destructor.
    int numPages; ///< Page count recorded by calculateNumPages().
    int Width, ///< Cached page width in pixels.
        Height, ///< Cached page height in pixels.
        samplesPerPixel, ///< Cached number of samples per pixel.
        bitsPerSample; ///< Cached bit depth of each sample.
    Constraints constraints; ///< Normalized crop bounds used by readinfo().

};
