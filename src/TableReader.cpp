//
// Created by dtaliun on 1/31/20.
//

#include "TableReader.h"

const unsigned int TableReader::DEFAULT_BUFFER_SIZE = 16777216u;

TableReader::TableReader(unsigned int buffer_size) noexcept(false): buffer(nullptr), max_line_length(0u), gzfile(nullptr) {
//    if (strcmp(zlibVersion(), ZLIB_VERSION) != 0) {
//        throw runtime_error("Incompatible ZLIB version");
//    }
    try {
        buffer = new char[buffer_size];
        max_line_length = buffer_size - 1u;
    } catch (bad_alloc& e) {
        throw runtime_error("Error in memory allocation");
    }
}

TableReader::~TableReader() {
    delete[] buffer;
    buffer = nullptr;
    if (gzfile != nullptr) {
        gzclose(gzfile);
        gzfile = nullptr;
    }
}

void TableReader::set_file_name(const string& file_name) {
    if (!is_open()) {
        this->file_name = file_name;
    }
}

const string& TableReader::get_file_name() {
    return file_name;
}

void TableReader::open() noexcept(false) {
    if (gzfile == nullptr) {
        gzfile = gzopen(file_name.c_str(), "rb");
        if (gzfile == nullptr) {
            throw runtime_error("Error while opening '" + file_name + "' file.");
        }
    }
}

void TableReader::close() noexcept(false) {
    if (gzfile != nullptr) {
        int gzerrno = 0;
        gzerrno = gzclose(gzfile);
        gzfile = nullptr;
        if (gzerrno != Z_OK) {
            throw runtime_error("Error while closing '" + file_name + "' file.");
        }
    }
}

void TableReader::reset() noexcept(false) {
    if (gzseek(gzfile, 0L, SEEK_SET) < 0) {
        throw runtime_error("Error while resetting '" + file_name + "' file.");
    }
}

bool TableReader::eof() {
    return gzeof(gzfile) > 0;
}

bool TableReader::is_open() {
    return gzfile != nullptr;
}

void TableReader::get_dim(int &nrow, int &ncol, char separator) noexcept(false) {
    vector<pair<int, int>> rows(1, pair<int, int>(1, 0)); // for each row we will store number of detected columns and number of characters.
    int c = 0;
    nrow = 0;
    ncol = 0;

    reset(); // make sure we are at the beginning of the file.
    while (!eof()) {
        while ((c = gzgetc(gzfile)) >= 0) {
            if ((char)c == '\n') {
                rows.emplace_back(1, 0);
            } else if ((char)c == separator) {
                rows.back().first += 1;
            } else {
                rows.back().second += 1;
            }
        }
        if ((c < 0) && (!eof())) {
            throw runtime_error("Error while reading '" + file_name + "' file.");
        }
    }
    for (auto&& row: rows) {
        if ((row.first == 1) && (row.second == 0)) { // empty row - we can ignore it
            continue;
        }
        if (row.first > row.second) { // some columns are empty
            throw runtime_error("Empty columns in '" + file_name + "' file.");
        }
        if (nrow == 0) {
            ncol = row.first;
        } else if (ncol != row.first) { // number of columns doesn't match between rows;
            throw runtime_error("Number of columns is different across rows in '" + file_name + "' file.");
        }
        ++nrow;
    }
    reset(); // move to the beginning of the file.
}

long int TableReader::read_row(vector<string>& tokens, char separator) noexcept(false) {
    long int i = 0;
    int c = 0;
    const char* token = buffer;

    tokens.clear();
    while ((i < max_line_length) && ((c = gzgetc(gzfile)) >= 0)) {
        buffer[i] = (char)c;

        if (buffer[i] == separator) {
            buffer[i] = '\0';
            tokens.emplace_back(token);
            token = buffer + i + 1;
        } else if (buffer[i] == '\n') {
            buffer[i] = '\0';
            tokens.emplace_back(token);
            return i;
        } else if (buffer[i] == '\r') {
            buffer[i] = '\0';
            tokens.emplace_back(token);
            if ((c = gzgetc(gzfile)) >= 0) {
                if ((char)c != '\n') {
                    c = gzungetc(c, gzfile);
                }
            }
            return i;
        }
        i += 1;
    }

    if ((c < 0) && (gzeof(gzfile) < 1)) {
        throw runtime_error("Error while reading rows from '" + file_name + "' file.");
    }

    return (i == 0 ? -1 : i);
}

int TableReader::check_format(int header_rows, int index_columns, int required_data_rows, int required_data_cols,  Format format, string& message) noexcept(false) {
    this->reset();
    stringstream fmt_message;
    vector<string> tokens;
    int nrow = 0;
    size_t pos = 0;
    float value = 0;

    while (this->read_row(tokens, '\t') >= 0) {
        ++nrow;
        if (nrow <= header_rows) { // Skip non-data rows
            continue;
        }
        if (tokens.size() != required_data_cols + index_columns) {
            fmt_message << "Error: incorrect number of columns in row " << nrow << " in '" << this->file_name << "'." << endl;
            message = fmt_message.str();
            return 0; // incorrect number of columns or empty row
        }
        for (unsigned int i = index_columns; i < tokens.size(); ++i) {
            if (tokens[i].length() == 0) {
                fmt_message << "Error: empty field at (row " << nrow << ", column " << i << ") in the file '" << this->file_name << "'." << endl;
                message = fmt_message.str();
                return 0;
            } else {
                switch (format) {
                    case Format::DIPLOID_GT:
                        if (tokens[i].compare("-9") != 0) {
                            if ((tokens[i].compare("0") != 0) && (tokens[i].compare("1") != 0) && (tokens[i].compare("2") != 0)) {
                                fmt_message << "Error: invalid value '" << tokens[i] << "' in (row " << nrow << ", column " << i << ") in the file '" << this->file_name << "'." << endl;
                                message = fmt_message.str();
                                return 0; // incorrect column value;
                            }
                        }
                        break;
                    case Format::NPLOID_GT:
                        if (tokens[i].compare("-9") != 0) {
                            if (tokens[i].find_first_not_of("0123456789") != string::npos) {
                                fmt_message << "Error: invalid value '" << tokens[i] << "' in (row " << nrow << ", column " << i << ") in the file '" << this->file_name << "'." << endl;
                                message = fmt_message.str();
                                return 0; // incorrect column value;
                            }
                        }
                        break;
                    case Format::FLOAT:
                        try {
                            value = stof(tokens[i], &pos);
                            if ((pos != tokens[i].length()) || (isinf(value))) {
                                throw runtime_error("");
                            }
                        } catch (...) {
                            fmt_message << "Error: invalid value '" << tokens[i] << "' in (row " << nrow << ", column " << i << ") in the file '" << this->file_name << "'." << endl;
                            message = fmt_message.str();
                            return 0;
                        }
                        break;
                }
            }
        }
    }

    if (nrow != (required_data_rows + header_rows)) {
        fmt_message << "Error: incorrect number of rows in the file '" << this->file_name << "'." << endl;
        message = fmt_message.str();
        return 0;
    }

    this->reset();
    message = "";
    return 1;
}