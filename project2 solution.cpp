
//I am making these changes over github for Lab #10

#include<string>
using std::string;
#include<vector>
using std::vector;
#include "bio.h"
//vectorization is function that recieve a string and return
//a vector<string> cuts in every 3 character, and ignore remainders
vector<string> vectorization(string sequence) {
    
    int length = sequence.size();
    
    vector<string> out_result(length/3);
    for (int x = 0; x < (length / 3); x++ ) {
    out_result[x] = sequence.substr(x*3, 3);
    }
    return out_result; 
}
//reverse the string order, for example, 
//if string is "asdf", it returns "fdsa" 
string reverse_sequence(string sequence) {
    int length = sequence.size();
    string reverse;
    for (int x = 0; x < length; x++) {
        reverse += sequence.back();
		sequence.pop_back();
    }
    return reverse;
}
//finding complement of the sequence, second parameter defines whether if
//is RNA complement or DNA complement
string complement_sequence(string sequence, string DNA_RNA) {
    string output;
    int length = sequence.size();
    if (DNA_RNA == "DNA") {
        for (int x = 0; x < length; x++) {
            if (sequence[x] == 'A') {
                output += 'T';
            }
            else if (sequence[x] == 'C') {
                output += 'G';
            }
            else if (sequence[x] == 'T') {
                output += 'A';
            }
            else if (sequence[x] == 'G') {
                output += 'C';
            }
        }
       
    }
    if (DNA_RNA == "RNA") {

        for (int x = 0; x < length; x++) {
            if (sequence[x] == 'A') {
                output += 'U';
            }
            if (sequence[x] == 'C') {
                output += 'G';
            }
            // in case we are finding complement of RNA sequence
            if (sequence[x] == 'U') {
                output += 'A';
            }
            if (sequence[x] == 'T') {
                output += 'A';
            }
            if (sequence[x] == 'G') {
                output += 'C';
            }
        }
    
    }
    return output;
}
// check whether it is valid sequence
bool IsValidDNASequence(const string& input) {
    int length = input.size();
    for (int x = 0; x <length; x++) {
        if (input[x] == 'A') {
            continue;
        }
        else if (input[x] == 'C') {
            continue;
        }
        else if (input[x] == 'T') {
            continue;
        }
        else if (input[x] == 'G') {
            continue;
        }
        else {
            return false;
        }
    }
    return true;
}


//getting reverse complement sequence
void GetReverseComplementSequence(const string & input, string* const output) {
    string reversed = reverse_sequence(input);
    *output = complement_sequence(reversed, "DNA");
}

// getting RNA transcript
string GetRNATranscript(const string & input){
    string output;
    string reversed = reverse_sequence(input);
    output = complement_sequence(reversed, "RNA");
    return output;
}

//getting 6 framses of codon
vector<vector<string>> GetReadingFramesAsCodons(const string& input) {
    vector<vector<string>> out_result(6);
    string DNA_copy = input;
    
    string RNA_copy = GetRNATranscript(DNA_copy);
    string original_1 = RNA_copy.substr(1);
    string original_2 = RNA_copy.substr(2);
    string anti_par_0 = (complement_sequence(reverse_sequence(RNA_copy), "RNA"));
    string anti_par_1 = anti_par_0.substr(1);
    string anti_par_2 = anti_par_0.substr(2);
    
    out_result[0] = (vectorization(RNA_copy));
    out_result[1] = (vectorization(original_1));
    out_result[2] = (vectorization(original_2));
    out_result[3] = (vectorization(anti_par_0));
    out_result[4] = (vectorization(anti_par_1));
    out_result[5] = (vectorization(anti_par_2));
    return out_result;
}
// translation of RNA sequence to amino acids
string Translate(const vector<string>& codon_sequence) {
    string output_string;
    int length = codon_sequence.size();
    vector<string> dictionary_sequence = { "GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
"AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
"GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
"UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
"CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
"ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
"GUG", "UAG", "UGA", "UAA" };
    int size_dictionary = dictionary_sequence.size();
    vector<string> dictionary_meaning = { "A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D",
"C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
"I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
"P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
"Y", "V", "V", "V", "V", "*", "*", "*" };
    for (int x = 0; x < length; x++) {
        for (int y = 0; y < size_dictionary; y++) {
            if (codon_sequence[x] == dictionary_sequence[y]) {
                output_string += dictionary_meaning[y];
            }
        }
    }
    
    return output_string;
}
// get longest frame
string GetLongestOpenReadingFrame(const string& DNA_sequence) {
    string longest;
    string amino_0;
    string amino_1;
    bool turn_token = 0;
    string total;
    int length_total;
    int longest_size = 0;
    int amino_0_size;
    for (vector<string> frame : GetReadingFramesAsCodons(DNA_sequence)) {
        amino_0 = "";
        total = Translate(frame);
        length_total = total.size();
        for (int x = 0; x < length_total; x++) {
            longest_size = longest.size();
            if (total[x] == 'M') {
                turn_token = 1;
            }
            if (turn_token == 1) {
                amino_0 += total[x];
            }
            if (total[x] == '*') {
                amino_0_size = amino_0.size();
                if (longest_size < amino_0_size) {
                    longest = amino_0;
                }
                amino_0 = "";
                turn_token = 0;
            }
        }
    }
    return longest;
}


