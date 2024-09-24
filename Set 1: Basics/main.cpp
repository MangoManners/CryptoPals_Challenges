#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>

#include "convertions.h"

void printMenu();
int selectOption();

void readFromUser(std::string, std::string&);
std::string getInputFileName();
std::string getOutputFileName();
void readFromFile(std::string&, std::string);
void writeToFile(std::string&, std::string);

std::string XOR(std::string&, std::string&);
std::string repeatingXOR(std::string&, std::string&);
std::vector<std::string> findSingleByte(std::string&);
std::vector<std::pair<char, int>> sortMap(std::map<char, int>);
int getScore(std::string);

void subBytes(matrix&);
void invSubBytes(matrix&);
void shiftRowByOneByteForward(unsigned char (&)[4]);
void shiftRowByOneByteBackward(unsigned char (&)[4]);
void shiftRows(matrix&);
void invShiftRows(matrix&);
unsigned char multiplicateInGaluaField(unsigned char, unsigned char);
void mixColumns(matrix&);
void invMixColumns(matrix&);

void subWord(word&);
void rotWord(word&);
void rCon(word&, int);
void XORWords(word&, word&);

void addRoundKey(matrix&, std::vector<word>&, int);

void expandKey(bytes&, std::vector<word>&, int);

void encryptBlock(bytes&, bytes&, std::vector<word>&, int);


void solveTask1();
void solveTask2();
void solveTask3();
void solveTask4();
void solveTask5();
void solveTask6();
void solveTask7();
void solveTask8();

int main() {
    printMenu();
    int opt = -1;
    while (opt != 0) {
        opt = selectOption();
        if (opt == 1) {
            solveTask1();
        } else if (opt == 2) {
            solveTask2();
        } else if (opt == 3) {
            solveTask3();
        } else if (opt == 4) {
            solveTask4();
        } else if (opt == 5) {
            solveTask5();
        } else if (opt == 6) {
            solveTask6();
        } else if (opt == 7) {
            solveTask7();
        } else if (opt == 8) {
            solveTask8();
        }
    }

    return 0;
}


void printMenu()
{
    std::cout << "Set 1: Basics. \n\t"
                "1. Convert hex to base64\n\t"
                "2. Fixed XOR\n\t"
                "3. Single-byte XOR cipher\n\t"
                "4. Detect single-character XOR\n\t"
                "5. Implement repeating-key XOR\n\t"
                "6. Break repeating-key XOR\n\t"
                "7. AES in ECB mode\n\t"
                "8. Detect AES in ECB mode\n\t"
                "0 - to exit\n";
}


int selectOption()
{
    int opt = -1;
    const std::vector<std::string> options {"0", "1", "2", "3", "4", "5", "6", "7", "8"};
    while (opt == -1) {
        std::string option = "";
        readFromUser("To select option enter its number: ", option);
        auto it = find(options.begin(), options.end(), option);
        if (it != options.end()) {
            opt = (it - options.begin());
        } else {
            std::cout << "Try again!\n";
        }
    }  

    return opt;
}


void readFromUser(std::string textToPrint, std::string& textToRead)
{
    std::cout << textToPrint;
    getline(std::cin, textToRead);
    //std::cout << textToRead;
}


std::string getInputFileName()
{
    std::string fileName = "";
    readFromUser("Enter file name to read from: ", fileName);
    if (fileName == "") {
        fileName = "input.txt";
    }
    return fileName;
}


std::string getOutputFileName()
{
    std::string fileName = "";
    readFromUser("Enter file name to write to: ", fileName);
    if (fileName == "") {
        fileName = "output.txt";
    }
    return fileName;
}


void readFromFile(std::string& text, std::string fileName)
{
    std::ifstream file (fileName);
    if (file.is_open()) {
        std::string temp = "";
        while (std::getline(file, temp)) {
            text += temp;
        }
        file.close();
    } else {
        std::cout << "Error during opening the file!\n"; 
    }
}


void writeToFile(std::string& text, std::string fileName)
{
    std::ofstream file (fileName, std::ios::out | std::ios::trunc);
    if (file.is_open()) {
        file << text;
        file.close();
    } else {
        std::cout << "Error during opening the file!\n"; 
    }
}



void solveTask1()
{
    std::string hex = "49276d206b696c6c696e6720796f757220627261696e206c696b65206120706f69736f6e6f7573206d757368726f6f6d";
    
    std::string encodedBase64 = "";
    convertHexToBase64(hex, encodedBase64);

    std::cout << "Result: " << encodedBase64 << "\n";
}


std::string XOR(std::string& a, std::string& b)
{
    while (a.length() < b.length()) {
        a = "0" + a;
    }
    while (a.length() > b.length()) {
        b = "0" + b;
    }
    bytes aBytes, bBytes;
    convertHexToBytes(a, aBytes);
    convertHexToBytes(b, bBytes);

    for (int i = 0; i < aBytes.size(); ++i) {
        aBytes[i] ^= bBytes[i];
    }

    a = "";
    convertBytesToHex(aBytes, a);

    return a;
}


void solveTask2()
{
    std::string a = "1c0111001f010100061a024b53535009181c";
    std::string b = "686974207468652062756c6c277320657965";
    std::string res = XOR(a, b);
    std::cout << "Result: " << res << "\n";
}


std::vector<std::pair<char, int>> sortMap(std::map<char, int> mp)
{
    std::vector<std::pair<char, int>> pairs; 
  
    for (std::pair<char, int> it : mp) { 
        pairs.push_back(it); 
    } 
  
    sort(pairs.begin(), pairs.end(), [](auto& a, auto& b) { return a.second > b.second; }); 

    return pairs;
}


int getScore(std::string ASCIIstr)
{
    std::vector<char> first13 {'e', 't', 'a', 'o', 'i', 'n', ' ', 's', 'h', 'r', 'd', 'l', 'u'};

    std::transform(ASCIIstr.begin(), ASCIIstr.end(), ASCIIstr.begin(), [](unsigned char c){ return std::tolower(c); });

    std::map<char, int> frequency;
    for (int i = 0; i < ASCIIstr.length(); ++i) { 
        if (frequency[ASCIIstr[i]]) {
            frequency[ASCIIstr[i]] += 1;
        } else {
            frequency[ASCIIstr[i]] = 1;
        }
    }

    std::vector<std::pair<char, int>> pairs = sortMap(frequency);
    std::vector<char> chars (13, '0');
    for (int i = 0; i < 13; ++i) {
        chars[i] = pairs[i].first;
    }

    int score = 0;
    for (char ch : chars) {
        int cnt = std::count(first13.begin(), first13.end(), ch); 
        if (cnt > 0) {
            score += frequency[ch];
        }
    }
    return score;
}


// {byte, string, score}
std::vector<std::string> findSingleByte(std::string& hex)
{
    int sizeRequired = hex.length();

    std::vector<std::string> bestResults (3, "0");
    int bestScore = 0;

    for (int i = 0; i < 256; ++i) {
        std::string singleByte = formatHex((byte)i);
        std::string byteStr = singleByte;
        
        while (byteStr.length() < sizeRequired) {
            byteStr += singleByte;
        }

        std::string xored = XOR(byteStr, hex);
        std::string result = "";
        convertHexToASCII(xored, result);

        if (isASCII(result)) {
            int score = getScore(result);
            if (score > bestScore) {
                bestResults[0] = singleByte;
                bestResults[1] = result;
                bestResults[2] = std::to_string(score);
                bestScore = score;
            }
        }           

    }
    return bestResults;
}


void solveTask3()
{
    std::string hex = "1b37373331363f78151b7f2b783431333d78397828372d363c78373e783a393b3736";
    std::vector<std::string> results = findSingleByte(hex);
    std::cout << "Result.\n\t" << 
                 "Byte: " << results[0] << "\n\t"
                 "String: " << results[1] << "\n\t"
                 "Score: " << results[2] << "\n";
}


void solveTask4()
{
    std::ifstream file (getInputFileName());
    std::string temp = "";
    if (!file.is_open()) { 
        std::cout << "Error during opening the file!\n"; 
        return;
    }
    std::vector<std::vector<std::string>> bestResults;
    while (std::getline(file, temp)) {
        std::vector<std::string> bestResult = findSingleByte(temp);
        bestResults.push_back(bestResult);
    } 
    file.close();

    int bestScore = 0;
    std::vector<std::string> bestRes;
    int score = 0;
    for (std::vector<std::string> res : bestResults) {
        score = std::stoi(res[2]);
        if (score > bestScore) {
            bestRes = res;
            bestScore = score;
        }
    }
    std::cout << "Result.\n\t" << 
                 "Byte: " << bestRes[0] << "\n\t"
                 "String: " << bestRes[1] << "\n\t"
                 "Score: " << bestRes[2] << "\n";
}


std::string repeatingXOR(std::string& plaintext, std::string& key)
{
    std::string ciphertext = "";
    int mod = key.length();
    int diff = plaintext.length() - mod;
    if (diff > 0) {
        for (int i = 0; i < diff; ++i) {
            key += key[i % mod];
        }
    }
    std::string plaintextHex = "";
    std::string keyHex = "";
    convertASCIIToHex(plaintext, plaintextHex);
    convertASCIIToHex(key, keyHex);
    ciphertext = XOR(plaintextHex, keyHex);

    return ciphertext;
}


void solveTask5()
{
    std::string text = "";
    readFromFile(text, getInputFileName());
    std::string key = "";
    readFromUser("Enter key: ", key);

    std::string xored = repeatingXOR(text, key);
    std::cout << "Result: " << xored << "\n";

}



void insertZeroes(std::string& str, int count)
{
    while (count != 0) {
        str = "0" + str;
        count -= 1;
    }

}


bytes slice(bytes& bytesVec, int start, int count)
{
    bytes slicedVec;
    while (count != 0) {
        slicedVec.push_back(bytesVec[start]);
        ++start;
        --count;
    }
    return slicedVec;
}


int computeHammingDistance(std::string binStr1, std::string binStr2)
{
    int dist = 0;

    int diff = binStr1.length() - binStr2.length();
    if (diff > 0) {
        insertZeroes(binStr2, diff);
    } else if (diff < 0) {
        insertZeroes(binStr1, diff*(-1));
    }

    int size = (binStr1.length() + binStr2.length()) / 2;

    for (int i = 0; i < size; ++i) {
        if (binStr1[i] != binStr2[i]) {
            dist += 1;
        }
    }

    return dist;
}


int guessKeyLength(int minLength, int maxLength, bytes& bytesVec)
{
    std::map<int, double> keyLengthProbabilities;
    std::vector<bytes> blocks (4);
    
    for (int keyLength = minLength; keyLength <= maxLength; ++keyLength) {
        blocks[0] = slice(bytesVec, 0, keyLength);
        blocks[1] = slice(bytesVec, keyLength, keyLength);
        blocks[2] = slice(bytesVec, keyLength*2, keyLength);
        blocks[3] = slice(bytesVec, keyLength*3, keyLength);
        double avgEditDist = 0.0;
        for (int i = 0; i < 4; ++i) {
            for (int j = i+1; j < 4; ++j) {
                std::string iBinary = "";
                std::string jBinary = "";
                convertBytesToBinary(blocks[i], iBinary);
                convertBytesToBinary(blocks[j], jBinary);
                avgEditDist += computeHammingDistance(iBinary, jBinary);
            }
        }
        avgEditDist /= keyLength;
        keyLengthProbabilities[keyLength] = avgEditDist;
    }
    std::vector<std::pair<int, double>> pairs; 
  
    for (std::pair<int, double> it : keyLengthProbabilities) { 
        pairs.push_back(it); 
    } 
  
    sort(pairs.begin(), pairs.end(), [](auto& a, auto& b) { return a.second < b.second; }); 
    
    return pairs[0].first;   
}


void splitIntoBlocks(bytes& orig, int size, std::vector<bytes>& blocks)
{
    int count = orig.size() / size;
    std::cout << orig.size() << " " << size << " " << count << "\n";
    for (int i = 0; i < count; ++i) {
        blocks.push_back(slice(orig, i * size, size));
    }
    int diff = orig.size() - (size * count);
    if (diff != 0) {
        bytes last (size, 0x00);
        int indLast = last.size() - 1;
        int indOrig = orig.size() - 1;
        while (diff > 0) {
            last[indLast] = orig[indOrig];
        }
    }

}


void transposeBlocks(std::vector<bytes>& orig, std::vector<bytes>& transposed)
{
    int rows = orig[0].size();
    int cols = orig.size();
    bytes temp (cols, 0x00);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            temp[col] = orig[col][row];
        }
        transposed.push_back(temp);
    }
}


void solveTask6()
{
    std::string text = "";
    std::string textASCII = "";
    readFromFile(text, getInputFileName());

    bytes ciphertextBytes; 
    convertBase64ToBytes(text, ciphertextBytes);
    convertBytesToASCII(ciphertextBytes, textASCII);

    int keySize = guessKeyLength(2, 40, ciphertextBytes);
    std::cout << "Guessed key length: " << keySize << "\n";

    std::vector<bytes> ciphertextBlocks;
    splitIntoBlocks(ciphertextBytes, keySize, ciphertextBlocks);

    std::vector<bytes> ciphertextTransposedBlocks;
    transposeBlocks(ciphertextBlocks, ciphertextTransposedBlocks);

    std::string key = "";
    std::string keyHex = "";
    for (bytes block : ciphertextTransposedBlocks) {
        std::string hexBlock = "";
        convertBytesToHex(block, hexBlock);
        keyHex += findSingleByte(hexBlock)[0];
    }

    convertHexToASCII(keyHex, key);
    std::cout << "Key: " << key << "\n";

    bytes keyBytes;
    convertASCIIToBytes(key, keyBytes);
    std::string plaintext = "";
    for (bytes& block : ciphertextBlocks) {
        for (int i = 0; i < block.size(); ++i) {
            block[i] ^= keyBytes[i];
        }
        for (unsigned char ch : block) {
            plaintext += (char)ch;
        }
    }
    writeToFile(plaintext, getOutputFileName());
}


void subBytes(matrix& state)
{
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            state[i][j] = S_BOX[(int) state[i][j]];
        }
    }
}


void invSubBytes(matrix& state)
{
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            state[i][j] = (unsigned char) (find(S_BOX.begin(), S_BOX.end(), state[i][j]) - S_BOX.begin());
        }
    }
}


void shiftRowByOneByteForward(unsigned char (&row)[4])
{
    for (int i = 0; i < 3; ++i) {
        unsigned char temp = row[i + 1];
        row[i + 1] = row[i];
        row[i] = temp;
    }

}

void shiftRowByOneByteBackward(unsigned char (&row)[4])
{
    for (int i = 3; i > 0; --i) {
        unsigned char temp = row[i - 1];
        row[i - 1] = row[i];
        row[i] = temp;
    }

}


void shiftRows(matrix& state)
{
    for (int i = 1; i < 4; ++i) {
        for (int count = i; count > 0; --count) {
            shiftRowByOneByteForward(state[i]);
        }
    }
}


void invShiftRows(matrix& state)
{
    for (int i = 1; i < 4; ++i) {
        for (int count = i; count > 0; --count) {
            shiftRowByOneByteBackward(state[i]);
        }
    }
}


void mixColumns(matrix& state)
{
    matrix tempState;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            tempState[i][j] = state[i][j];
        }
    }
    for (int i = 0; i < 4; ++i) { 
        for (int j = 0; j < 4; ++j) {
            unsigned char temp = 0x00;
            //std::cout << i << ", " << j << ": \n";
            for (int row = 0; row < 4; ++row) {
                //std::cout << formatHex(MIX[i][row]) << " * " << formatHex(tempState[row][j]);
                temp ^= multiplicateInGaluaField(MIX[i][row], tempState[row][j]);        
                //std::cout << " = " << formatHex(temp) << "\n";
            }
            state[i][j] = temp;
        }
    }
}


void invMixColumns(matrix& state)
{
    matrix tempState;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            tempState[i][j] = state[i][j];
        }
    }
    for (int i = 0; i < 4; ++i) { 
        for (int j = 0; j < 4; ++j) {
            unsigned char temp = 0x00;
            //std::cout << i << ", " << j << ": \n";
            for (int row = 0; row < 4; ++row) {
                //std::cout << formatHex(MIX[i][row]) << " * " << formatHex(tempState[row][j]);
                temp ^= multiplicateInGaluaField(INV_MIX[i][row], tempState[row][j]);        
                //std::cout << " = " << formatHex(temp) << "\n";
            }
            state[i][j] = temp;
        }
    }
}


unsigned char multiplicateInGaluaField(unsigned char a, unsigned char b)
{
    unsigned char c = 0x00;
    while (b > 0) {
        if (b & 1) {
            c ^= a;
        }
        a = (a << 1) ^ (a & 0x80 ? 0x1b : 0x00);
        b >>= 1;
    }
    return c;
}


void subWord(word& w)
{
    for (int i = 0; i < 4; ++i) {
        w[i] = S_BOX[w[i]];
    }
}


void rotWord(word& w)
{
    shiftRowByOneByteForward(w);
}


void rCon(word& w, int ind)
{
    unsigned char firstByte = 0x01;
    for (int i = 1; i < ind; ++i) {
        firstByte = multiplicateInGaluaField(firstByte, 0x02);
    }
    w[0] = firstByte;
    for (int i = 1; i < 4; ++i) {
        w[i] = 0x00;
    }
}


void XORWords(word& a, word& b)
{
    for (int i = 0; i < 4; ++i) {
        a[i] ^= b[i];
    }
}


void copyWord(word& from, word& to)
{
    for (int i = 0; i < 4; ++i) {
        to[i] = from[i];
    }
}


void expandKey(bytes& key, std::vector<word>& keySchedule, int Nk)
{
    for (int i = 0; i < Nk; ++i) {
        for (int j = 0; j < 4; ++j) {
            keySchedule[i][j] = key[4 * i + j];
        }
    }
    
    word r {0x00, 0x00, 0x00, 0x00};
    for (int i = Nk; i < keySchedule.size(); ++i) {
        //std::cout << "Iter: " << i << "\n";
        word temp {0x00, 0x00, 0x00, 0x00};
        copyWord(keySchedule[i - 1], temp);
        //printWord(temp);
        if (i % Nk == 0) {
            rotWord(temp);
            //printWord(temp);
            subWord(temp);
            //printWord(temp);
            rCon(r, i / Nk);
            //printWord(r);
            XORWords(temp, r);
            //printWord(temp);
        } else if (Nk > 6 && i % Nk == 4) {
            subWord(temp);
        }
        XORWords(temp, keySchedule[i - Nk]);
        //printWord(temp);
        copyWord(temp, keySchedule[i]);

    }
}


void addRoundKey(matrix& state, std::vector<word>& keys, int round) 
{
    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            state[row][col] ^= keys[round * 4 + col][row];
        }
    }
}



void encryptBlock(block& plaintext, block& ciphertext, std::vector<word>& keys, int rounds)
{
    matrix state;
    convertBlockToMatrix(plaintext, state);
    addRoundKey(state, keys, 0);

    for (int round = 1; round < rounds; ++round) {
        subBytes(state);
        shiftRows(state);
        mixColumns(state);
        addRoundKey(state, keys, round);
    }

    subBytes(state);
    shiftRows(state);
    addRoundKey(state, keys, rounds);

    convertMatrixToBlock(state, ciphertext);

}


void decryptBlock(block& ciphertext, block& plaintext, std::vector<word>& keys, int rounds)
{
    matrix state;
    convertBlockToMatrix(ciphertext, state);
    addRoundKey(state, keys, rounds);

    for (int round = 9; round > 0; --round) {
        invShiftRows(state);
        invSubBytes(state);
        addRoundKey(state, keys, round);
        invMixColumns(state);
    }

    invShiftRows(state);
    invSubBytes(state);
    addRoundKey(state, keys, 0);

    convertMatrixToBlock(state, plaintext);

}


void askParams(params& algorithm)
{
    std::string info = "Choose algorithm. \n\t1. AES-128 \n\t2. AES-192 \n\t3. AES-256\n"
                    "Enter num: ";
    int choice = 1;
    std::string opt = "";
    readFromUser(info, opt);
    
    if (opt == "2") {
        choice = 2;
    } else if (opt == "3") {
        choice = 3;
    }
    auto iter = encryptionAlgorithmParams.find(choice);
    if (iter != encryptionAlgorithmParams.end()) {
        algorithm.Nk = iter->second.Nk;
        algorithm.Nb = iter->second.Nb;
        algorithm.Nr = iter->second.Nr;
    }
    std::cout << "Next parameters will be used: Nk = " << algorithm.Nk << 
                                              ", Nb = " << algorithm.Nb <<
                                              ", Nr = " << algorithm.Nr << "\n";
}


void askKey(bytes& text)
{
    std::string keyStr = "YELLOW SUBMARINE";
    std::string keyHex = "";
    convertASCIIToHex(keyStr, keyHex);
    convertHexToBytes(keyHex, text);
}


void askMode(int& mode)
{
    std::string info = "Choose mode. \n\t0. encrypt \n\t1. decrypt\n"
                    "Enter num: ";
    std::string opt = "";
    readFromUser(info, opt);
    if (opt == "0") {
        mode = 0;
    } else if (opt == "1") {
        mode = 1;
    }
}


//mode = encrypt/decrypt
void setup(params& algorithm, bytes& text, bytes& key, int& mode)
{
    askParams(algorithm);
    std::string textBase64 = "";

    readFromFile(textBase64, getInputFileName());
    
    convertBase64ToBytes(textBase64, text);
    //printBytes(text);

    askKey(key);
    std::cout << "Your key in bytes.\n";
    printBytes(key);

    askMode(mode);
    std::cout << "You selected mode: " << mode << "\n";

}


void solveTask7()
{
    params algorithm;

    bytes textBytes;
    bytes keyBytes;

    std::string text = "";
    int mode = 0;
    setup(algorithm, textBytes, keyBytes, mode);

    std::vector<word> keySchedule (algorithm.Nb*(algorithm.Nr + 1));
    expandKey(keyBytes, keySchedule, algorithm.Nk); 

    for (int i = 0; i < textBytes.size(); i += 16) {
        std::string temp = "";
        bytes tempBytes;
        bytes blockIn = slice(textBytes, i, 16);
        bytes blockOut (blockIn.size());
        block aesBlockIn;
        block aesBlockOut;
        convertBytesToBlock(blockIn, aesBlockIn);
        convertBytesToBlock(blockOut, aesBlockOut);
        decryptBlock(aesBlockIn, aesBlockOut, keySchedule, algorithm.Nr);
        
        convertBlockToBytes(aesBlockOut, tempBytes);

        convertBytesToHex(tempBytes, temp);
        text += temp;
    }    
    std::string textASCII = "";
    convertHexToASCII(text, textASCII);
    writeToFile(textASCII, getOutputFileName());
}


bool compareBytes(const bytes& a, const bytes& b)
{
    for (int i = 0; i < 16; ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }

    return true;
}


void solveTask8()
{
    std::string fileName = getInputFileName();
    std::cout << fileName << "\n";
    std::ifstream file (fileName);
    std::vector<bytes> inputs;
    if (file.is_open()) {
        std::string temp = "";
        while (std::getline(file, temp)) {
            bytes tempBytes;
            convertHexToBytes(temp, tempBytes);
            inputs.push_back(tempBytes);
        }
        file.close();
    } else {
        std::cout << "Error during opening the file!\n"; 
    }

    std::cout << inputs.size() << "\n";

    std::vector<bytes> suitable;
    for (bytes input : inputs) {
        std::vector<bytes> blocks;
        bool flag = false;
        
        for (int i = 0; i < input.size(); i += 16) {
            bytes temp = slice(input, i, 16);
            
            for (int j = 0; j < blocks.size(); ++j) {
                if (compareBytes(temp, blocks[j]) && !flag) {
                    suitable.push_back(input);
                    flag = true;
                }
            }

            blocks.push_back(temp);
        }
    }

    for (int i = 0; i < suitable.size(); ++i) {
        printBytes(suitable[i]);
    }

}