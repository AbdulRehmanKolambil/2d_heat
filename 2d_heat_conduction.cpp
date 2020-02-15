#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;

int main()
{

    int n = 201;
    float x[n],y[n],Tnew[n][n],Told[n][n],Tprime[n][n],error[n][n];
	int i,j,iteration;
	float dx,dy,dt,endtime,timenow,diffus,maxerror,tolerance,k1,k2,SOR;

	/*---------------------------------------------------------------------------------
	n            :Number of grid points along x and y
	x[n]         :X co-ordinates of grid pints
	y[n]         :Y co-ordinates of grid pints
	Tnew[n][n]   :Newly calculated Temperature at time = endtime
	Told[n][n]   :Temperature at time = endtime - time step(dt)
	Tprime[n][n] :Temperature at endtime, and at previous iteration of space derivative
	-----------------------------------------------------------------------------------*/

    dx = 1/float(n);        //Grid size along x
    dy = 1/float(n);        //Grid size along x
    dt = 0.1;               //Time step
    timenow = 0.0;          //Time variable
    endtime = 7200.0;       //End time
    diffus = 0001.6563;     //Thermal Diffusivity assuming Aluminium
    maxerror = 100;         //Error initialization
    tolerance = 1e-5;       //Convergence Criterion
    SOR = 1.5;              //Successive Over-Relaxation Factor

	// Generating a uniformly structured grid
	for (i=0;i<n;i++)
    {
        x [i] = float(i)/(n-1);
        y [i] = float(i)/(n-1);
	}

    //Initializing the temperature
    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {
                if (pow((x[i]-0.5),2)+pow((y[j]-0.5),2)<0.2)
                {
                    Tnew[i][j]=40;
                    Told[i][j]=40;
                    Tprime[i][j]=40;
                }
                else
                {
                    Tnew[i][j]=20;
                    Told[i][j]=20;
                    Tprime[i][j]=20;
                }
        }
    }

    k1 = diffus*dt/pow(dx,2);
    k2 = diffus*dt/pow(dy,2);

    iteration =0;   //Iteration count

    //Solution using FDM
  while(timenow<(endtime+dt)) //Time loop
   {
        while(maxerror>tolerance) //Iterative loop to check convergence of spacial derivatives
        {
            for (i=1;i<n-1;i++){
                for (j=1;j<n-1;j++){
                    Tnew[i][j] =  (Told[i][j] + k1*(Tnew[i-1][j]+Tprime[i+1][j]) + k2*(Tnew[i][j-1]+Tprime[i][j+1]))/(1+2*k1+2*k2);
                }
            }
            iteration++;

            //Calculating the maximum error
            for (i=1;i<n;i++)
            {
                for (j=1;j<n;j++)
                {
                    error[i][j] = abs(Tnew[i][j]-Tprime[i][j]);
                    if (error[i][j]>error[i-1][j-1]){maxerror = error[i][j];}
                    else {maxerror = error[i-1][j-1];}
                }
            }

            //Output message
            cout << "Time =" << timenow << "\t\tIteration = " << iteration << "\tError = " << maxerror << endl;

            //Updating Tprime after each iteration
            for (i=1;i<n-1;i++){
                for (j=1;j<n-1;j++){
                   Tprime[i][j] = SOR*Tnew[i][j] + (1-SOR)*Tprime[i][j];
                }
           }

        }
         //Updating Tprime after each time-step
        for (i=1;i<n-1;i++){
                for (j=1;j<n-1;j++){
                   Told[i][j] = Tnew[i][j];
                }
        }

        //Increasing time-step and initializing Error
        timenow = timenow + dt;
        maxerror = 100;
    }

    //Writing "x co-ordinate,y-co-ordinate,Temperature" at time = endtime to a CSV file
    ofstream myfile;
    myfile.open ("Temperature.csv");
    for (i=0;i<n;i++)
            {
                for (j=0;j<n;j++)
                {
                    myfile << x[i] << "," << y[j] << "," << Tnew[i][j] << "\n";
                }
            }
    myfile.close();

    return 0;
}
