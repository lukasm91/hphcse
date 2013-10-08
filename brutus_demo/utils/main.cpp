#include <iostream>
#include "ArgumentParser.h"
#include "Timer.h"

using namespace std;

int main (int argc, char *argv[])
{
	ArgumentParser parser(argc, argv);

	const int nSteps = parser("-steps").asInt(10);	// number of steps, default 10
	const string msg = parser("-msg").asString("Default message");	// message to print within the loop
	const double num = parser("-num").asDouble(0.0);	// just a floating point number
	const bool flag = parser("-print").asBool(1);		// 

	cout << "Floating point number = " << num << endl;		

	Timer timer; // declaration of the timer
		
	timer.start();	// start the timer
	for (int i = 0; i < nSteps; i++)
	{
		if (flag)	// check if print flag was set
		{
				cout << msg << ", iter:" << i << endl;	// print message and iteration number
		}
		sleep(1);	// sleep for 1 second
	}
	double elapsed_time = timer.stop();	// stop the timer and get the elapsed time
	
	cout << "Elapsed time: " << elapsed_time << " seconds\n";  	

	return 0;
}
