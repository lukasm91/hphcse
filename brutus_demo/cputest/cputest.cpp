/* Simple cpu code that prints the hostname of the current node */
#include <stdio.h>
#include <unistd.h>

int main(int  argc, char *argv[])
{
	char name[256];
	
	gethostname(name, 256);
	
	printf("Hello, I am running on host: %s\n", name);
	
	return 0;
}
