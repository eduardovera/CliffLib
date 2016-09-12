#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

#include <iostream>
#include <map>
#include <bitset>
#include <algorithm>

using namespace std;

class multivector {

    private:
        std::map<uint, double> M;

        string get_base(int index) {
            std::string binary = std::bitset<32>(index).to_string();
            reverse(binary.begin(), binary.end());
            string ret = "";
            string s = "1";

            size_t pos = binary.find(s);
            while(pos != string::npos) {
                ret.append("e");
                ret.append(to_string(pos + 1));
                pos = binary.find(s, pos+1);
                if (pos != string::npos) {
                    ret.append("^");
                }
            }
            return ret;
        }

    public:

        multivector(std::initializer_list<double> l) {
            int index = 0;
            for (std::initializer_list<double>::iterator it = l.begin(); it != l.end(); ++it) {
                if ((*it) != 0.0) {
                    M[index] = (*it);
                }
                index++;
            }
        }

        void print() {
            for (std::map<uint, double>::iterator it = M.begin(); it != M.end(); ++it) {
                cout << ((*it).second > 0.0 ? "+" : "") << (*it).second << "(" << get_base((*it).first) << ")";
            }
            cout << endl;
        }
};

#endif // MULTIVECTOR_H
