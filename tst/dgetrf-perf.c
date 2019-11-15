#include "dgetrf.h"
#include "algonum.h"

void test_result() {

}

int main(int argc, char **argv) {
	char *str[] = {"OK", "NOK"};

	int r1 = !!testall_dgetrf( my_dgetrf );

	test_result();
	return 0;
}
