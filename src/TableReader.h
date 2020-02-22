//
// Created by dtaliun on 1/31/20.
//

#ifndef LASER_TABLEREADER_H
#define LASER_TABLEREADER_H

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <cmath>
#include <zlib.h>

using namespace std;

class TableReader {

private:
    char* buffer;
    unsigned int max_line_length;

    string file_name;
    gzFile gzfile;

public:
    enum Format { DIPLOID_GT, NPLOID_GT, FLOAT, SEQ };

    static const unsigned int DEFAULT_BUFFER_SIZE;

    TableReader(unsigned int buffer_size = DEFAULT_BUFFER_SIZE) noexcept(false);
    virtual ~TableReader();

    void set_file_name(const string& file_name);
    const string& get_file_name();

    void open() noexcept(false);
    void close() noexcept(false);
    void reset() noexcept(false);
    bool eof();
    bool is_open();

    void get_dim(int &nrow, int &ncol, char separator) noexcept(false);
    long int read_row(vector<string>& tokens, char separator) noexcept(false);
    int check_format(int header_rows, int index_columns, int required_data_rows, int required_data_cols,  Format format, string& message) noexcept(false);
};


#endif //LASER_TABLEREADER_H
