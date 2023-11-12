#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <cstdlib>
#include <cstdio>
#include <iostream>

#define ASSERT(cond, msg) if(cond){} else { \
	char buff[128];							\
	sprintf(buff, "Assetion Failed [LINE %d, FILE %s]: \n\t%s", __LINE__, __FILE__, msg); 	\
	std::cerr << buff << std::endl;			\
	exit(1); } 								\

#define SHORT_ASSERT(cond) ASSERT(cond, "")

#define ERROR(msg) ASSERT(0 == 0, msg)

#define UNIT_TEST(cond, msg) if(cond){ \
		std::cout << "Test SUCCESSEDED: " << #cond << "\n\n" << std::endl;} \
		else {std::cout << "Test Failed [LINE " << __LINE__ << ", FILE " << __FILE__ << "]: \n\t" << msg << std::endl; exit(1);}


#define printvar(x) std::cout << #x << " = " << x << std::endl;

#endif // UTILS_H_INCLUDED
