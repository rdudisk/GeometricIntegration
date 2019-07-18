#include <iostream>

using namespace std;

class A 
{
	public:
		static const unsigned int DOF = 3;

};

template <typename A>
class TB
{
	private:
		static const int i = A::DOF;

	public:
		void print ( )
		{
			std::cout << "DOF = " << i << std::endl;
		}
};

typedef TB<A> B;

int main (int argc, char* argv[]) {
	B b;
	b.print();

	return 0;
}
