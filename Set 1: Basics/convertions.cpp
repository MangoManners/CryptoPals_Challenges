#include "convertions.h"

void padding(std::string& str, int mod, int side = 0)
{
    if (side == 0) {
        while (str.length() % mod != 0) {
            str = "0" + str;
        }
    } else {
        while (str.length() % mod != 0) {
            str += "0";
        }
    }
}


void paddingBytes(bytes& vec, int reqSize)
{
    while (vec.size() < reqSize) {
        vec.insert(vec.begin(), 0x00);
    }
}
    

std::string formatHex(byte ch)
{
    std::string hex = "";
    std::stringstream ss;
    ss << std::hex << (int)ch;
    hex += ss.str();
    padding(hex, 2);
    return hex;
}


void printBytes(bytes& arr)
{
    std::cout << "Bytes: \n\t";
    for (unsigned char byte : arr) {
        std::cout << formatHex(byte) << " ";
    }
    std::cout << "\n";
}


// Base64 <-> Bytes
void convertBase64ToBytes(std::string& encoded, bytes& raw)
{
    //remove =
    encoded.erase(std::remove(encoded.begin(), encoded.end(), '='), encoded.end());
    //std::cout << encoded << "\n";
    //Convert to binary string
    std::string binary = "";
    convertBase64ToBinary(encoded, binary);
    //Convert to bytes
    convertBinaryToBytes(binary, raw);
}


void convertBytesToBase64(bytes& raw, std::string& encoded)
{
    encoded = "";
    std::string binary = "";
    convertBytesToBinary(raw, binary);
    convertBinaryToBase64(binary, encoded);
}


// ASCII <-> Bytes
void convertASCIIToBytes(std::string& text, bytes& raw)
{
    raw.clear();
    for (char ch : text) {
        raw.push_back(ch);
    }
}


void convertBytesToASCII(bytes& raw, std::string& text)
{
    text = "";
    for (byte b : raw) {
        text += (char)b;
    }
}


// Hex <-> Bytes
void convertHexToBytes(std::string& str, bytes& vec)
{
    vec.clear();
    std::string binary = "";
    convertHexToBinary(str, binary);
    convertBinaryToBytes(binary, vec);
}


void convertBytesToHex(bytes& vec, std::string& str)
{
    str = "";
    
    for (byte b : vec) {
        str += formatHex(b);
    }
}


// Binary <-> Bytes
void convertBinaryToBytes(std::string& str, bytes& vec)
{
    vec.clear();
    padding(str, 8);

    for (int i = 0; i < str.length(); i += 8) {
        std::bitset<8> bits (str.substr(i, 8));
        vec.push_back((unsigned char) bits.to_ulong());
    }

}


void convertBytesToBinary(bytes& vec, std::string& str)
{
    str = "";
    
    for (byte b : vec) {
        std::bitset<8> bits (b);
        str += bits.to_string();
    }

}


// Base64 <-> Binary
void convertBase64ToBinary(std::string& encoded, std::string& binary)
{
    binary = "";
    for (char ch : encoded) {
        long unsigned ind = base64Chars.find(ch);
        std::bitset<6> section{ind};
        binary += section.to_string();
    }
    //delete last zeroes
    while (binary[binary.length() - 1] == '0') {
        binary.pop_back();
    }
}


void convertBinaryToBase64(std::string& binary, std::string& encoded)
{
    encoded = "";

    padding(binary, 6, 1);

    for (int i = 0; i < binary.length(); i += 6) {
        std::bitset<6> charInd {binary.substr(i, 6)};
        encoded += base64Chars[charInd.to_ulong()];
    }

    int blocks = binary.length() / 6;

    if (blocks % 2 == 0 && blocks % 4 != 0) {
        encoded += "==";
    } else if (blocks % 3 == 0) {
        encoded += "=";
    }
}


// Hex <-> Binary
void convertBinaryToHex(std::string& binary, std::string& hex)
{
    hex = "";
    padding(binary, 4);
    for (int i = 0; i < binary.length(); i += 4) {
        std::string bits = binary.substr(i, 4);
        std::bitset<4> hexChar (bits);
        std::stringstream ss;
        ss << std::hex << hexChar.to_ulong();
        hex += ss.str();
    }
}


void convertHexToBinary(std::string& hex, std::string& binary)
{
    binary = "";
    for (int i = 0; i < hex.length(); ++i) {
        std::bitset<4> b {(unsigned int)(0x00 ^ ((hex[i] > '9') ? hex[i] - 87 : hex[i]))};
        binary += b.to_string();
    }
}


// Base64 <-> Hex
void convertBase64ToHex(std::string& encoded, std::string& hex)
{
    hex = "";
    std::string binary = "";
    convertBase64ToBinary(encoded, binary);
    convertBinaryToHex(binary, hex);
}


void convertHexToBase64(std::string& hex, std::string& encoded)
{
    encoded = "";
    std::string binary = "";
    convertHexToBinary(hex, binary);
    convertBinaryToBase64(binary, encoded);
}


// ASCII <-> Hex
void convertHexToASCII(std::string& hex, std::string& ASCII)
{
    ASCII = ""; 
    std::string part = "";
    for (int i = 0; i < hex.length(); i += 2) {
        part = hex.substr(i, 2);
        char ch = std::stoul(part, nullptr, 16);
        ASCII += ch;
    }
}

void convertASCIIToHex(std::string& ASCII, std::string& hex)
{
    hex = "";
    std::string byte = "";
    for (char ch : ASCII) {
        hex += formatHex(ch);
    }
}


// Bytes <-> Block
void convertBytesToBlock(bytes& vec, block& arr)
{
    int size = vec.size();
    if (size > 16) {
        std::cout << "Wrong size.\n";
    }
    paddingBytes(vec, 16);
    for (int i = 0; i < 16; ++i) {
        arr[i] = vec[i];
    } 
}


void convertBlockToBytes(block& arr, bytes& vec)
{
    vec.clear();
    for (int i = 0; i < 16; ++i) {
        vec.push_back(arr[i]);
    }
}


// Block <-> Matrix
void convertBlockToMatrix(block& arr, matrix& state)
{
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            state[i][j] = arr[4 * j + i];
        }
    }
}


void convertMatrixToBlock(matrix& state, block& arr)
{
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            arr[4 * j + i] = state[i][j];
        }
    }
}


bool isHex(std::string& text)
{
    for (char ch : text) {
        bool flag = false;
        for (char r : hexChars) {
            if (ch == r) {
                flag = true;
            }
        }
        if (!flag) {
            return false;
        }
    }

    return true;
}

bool isBase64(std::string& text)
{
    for (char ch : text) {
        bool flag = false;
        for (char r : base64Chars) {
            if (ch == r) {
                flag = true;
            }
        }
        if (!flag) {
            return false;
        }
    }

    return true;
}


bool isASCII(std::string& text)
{
    for (int c : text) {
        if (c > 127 || c < 0 ) {
            return false;
        }
    }
    return true;
}

