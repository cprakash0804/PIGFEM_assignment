#include <fstream>
#include <vector>

using namespace std;



int main (void)
{
	vector<float> vec = {0.0, 1.0, 2.0, 3.0};
	ofstream out;
	out.open("data.dat", ofstream::out | ofstream::binary);
	out.write(reinterpret_cast<char *>(&vec[0]), sizeof(float)*vec.size());
	out.write(reinterpret_cast<char *>(&vec[0]), sizeof(float)*vec.size());
	out.close();
}
