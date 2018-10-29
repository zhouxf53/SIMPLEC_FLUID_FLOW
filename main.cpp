#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream> 
#include <string> 
#include <algorithm>
#include <iomanip>
using namespace std;
const int IPCV = 40;              //row
const int JPCV = 120;


void TDMAI(double xx[IPCV], double ax[IPCV], double bx[IPCV], double cx[IPCV], double dx[IPCV])
{

	double u[IPCV], y[IPCV];
	u[IPCV - 1] = 0;
	u[0] = cx[0] / bx[0];
	y[0] = dx[0] / bx[0];
	int i;
	for (i = 1; i < IPCV - 1; i++)
	{
		u[i] = cx[i] / (bx[i] - u[i - 1] * ax[i]);            //simplify the equation

	}
	for (i = 1; i < IPCV; i++)
	{
		y[i] = (dx[i] - y[i - 1] * ax[i]) / (bx[i] - u[i - 1] * ax[i]);
	}

	//solve it

	xx[IPCV - 1] = y[IPCV - 1];
	for (i = IPCV - 2; i >= 0; i--)
	{
		xx[i] = y[i] - u[i] * xx[i + 1];
	}


}
void TDMAI1(double x[IPCV - 1], double a[IPCV - 1], double b[IPCV - 1], double c[IPCV - 1], double d[IPCV - 1])
{
	//b[i]*Tp=-a[i]*Tw-c[i]*Te+d[i]
	//double a[IPCV-1],b[IPCV-1],c[IPCV-1],d[IPCV-1];
	double u[IPCV - 1], y[IPCV - 1];
	u[IPCV - 1 - 1] = 0;
	u[0] = c[0] / b[0];
	y[0] = d[0] / b[0];
	int i;
	for (i = 1; i < IPCV - 1 - 1; i++)
	{
		u[i] = c[i] / (b[i] - u[i - 1] * a[i]);            //simplify the equation

	}
	for (i = 1; i < IPCV - 1; i++)
	{
		y[i] = (d[i] - y[i - 1] * a[i]) / (b[i] - u[i - 1] * a[i]);
	}

	//solve it

	x[IPCV - 1 - 1] = y[IPCV - 1 - 1];
	for (i = IPCV - 1 - 2; i >= 0; i--)
	{
		x[i] = y[i] - u[i] * x[i + 1];
	}


}
void TDMAJ(double xx[JPCV], double ax[JPCV], double bx[JPCV], double cx[JPCV], double dx[JPCV])
{

	double u[JPCV], y[JPCV];
	u[JPCV - 1] = 0;
	u[0] = cx[0] / bx[0];
	y[0] = dx[0] / bx[0];
	int i;
	for (i = 1; i < JPCV - 1; i++)
	{
		u[i] = cx[i] / (bx[i] - u[i - 1] * ax[i]);            //simplify the equation

	}
	for (i = 1; i < JPCV; i++)
	{
		y[i] = (dx[i] - y[i - 1] * ax[i]) / (bx[i] - u[i - 1] * ax[i]);
	}

	//solve it

	xx[JPCV - 1] = y[JPCV - 1];
	for (i = JPCV - 2; i >= 0; i--)
	{
		xx[i] = y[i] - u[i] * xx[i + 1];
	}


}
void TDMAJ1(double x[JPCV - 1], double a[JPCV - 1], double b[JPCV - 1], double c[JPCV - 1], double d[JPCV - 1])
{
	//b[i]*Tp=-a[i]*Tw-c[i]*Te+d[i]
	//double a[JPCV-1],b[JPCV-1],c[JPCV-1],d[JPCV-1];
	double u[JPCV - 1], y[JPCV - 1];
	u[JPCV - 1 - 1] = 0;
	u[0] = c[0] / b[0];
	y[0] = d[0] / b[0];
	int i;
	for (i = 1; i < JPCV - 1 - 1; i++)
	{
		u[i] = c[i] / (b[i] - u[i - 1] * a[i]);            //simplify the equation

	}
	for (i = 1; i < JPCV - 1; i++)
	{
		y[i] = (d[i] - y[i - 1] * a[i]) / (b[i] - u[i - 1] * a[i]);
	}

	//solve it

	x[JPCV - 1 - 1] = y[JPCV - 1 - 1];
	for (i = JPCV - 1 - 2; i >= 0; i--)
	{
		x[i] = y[i] - u[i] * x[i + 1];
	}


}
double max(double a, double b)
{
	if (a >= b)
	{
		return a;
	}
	else;
	return b;
}

void main()
{
	//int N=N;
	//call properties
	double rho = 997.0;
	double K = 0.563;
	double mu = 0.000871;
	double cp = 4179;
	double Re = 200;

	//ICS and BCS
	double p[IPCV + 2][JPCV + 2] = { 0 };
	double pc[IPCV + 2][JPCV + 2] = { 0 };
	double T[IPCV + 2][JPCV + 2] = { 0 };
	double Told[IPCV + 2][JPCV + 2] = { 0 };
	double dT[IPCV][JPCV] = { 0 };

	double u[IPCV + 2][JPCV + 1] = { 0 };
	double uold[IPCV + 2][JPCV + 1] = { 0 };

	double du[IPCV + 2][JPCV + 1] = { 0 };
	double v[IPCV + 1][JPCV + 2] = { 0 };

	double vold[IPCV + 1][JPCV + 2] = { 0 };
	double dv[IPCV + 1][JPCV + 2] = { 0 };
	double Vinf = Re * mu / rho / 0.04;
	int i, j;
	int id;
	//	int iu,ju,iv,jv,ip,jp;
	int k;
	for (i = 0; i < IPCV + 2; i++)
	{
		u[i][0] = Vinf;
		T[i][0] = 27;
	}
	for (j = 0; j < JPCV + 2; j++)
	{
		T[0][j] = 100;
		T[IPCV + 1][j] = 100;
	}
	double w = 0.01; //////////////relax factor for u and v
	double wp = 0.01; /////////////relax factor for p
	double Rp = 0;  ////these are residue for p, u, and v
	double Rv = 0;
	double Ru = 0;
	double temp = 0;
	double delx = 2.0 / JPCV;
	double dely = 0.02 / IPCV;



	//u coefficient
	double fw[IPCV][JPCV - 1], fe[IPCV][JPCV - 1], fs[IPCV][JPCV - 1], fn[IPCV][JPCV - 1];
	double ftw[IPCV][JPCV], fte[IPCV][JPCV], fts[IPCV][JPCV], ftn[IPCV][JPCV];  //temperature
	double dw[IPCV][JPCV - 1], de[IPCV][JPCV - 1], ds[IPCV][JPCV - 1], dn[IPCV][JPCV - 1];
	double Ps[IPCV][JPCV - 1], Pn[IPCV][JPCV - 1], Pw[IPCV][JPCV - 1], Pe[IPCV][JPCV - 1];
	double Pts[IPCV][JPCV], Ptn[IPCV][JPCV], Ptw[IPCV][JPCV], Pte[IPCV][JPCV];    //temperature
	double APs[IPCV][JPCV - 1], APn[IPCV][JPCV - 1], APw[IPCV][JPCV - 1], APe[IPCV][JPCV - 1];
	double APts[IPCV][JPCV], APtn[IPCV][JPCV], APtw[IPCV][JPCV], APte[IPCV][JPCV];  //temperature
	double a[JPCV - 1], b[JPCV - 1], c[JPCV - 1], d[JPCV - 1], d1[JPCV - 1], d2[JPCV - 1];
	double ax[IPCV], bx[IPCV], cx[IPCV], dx[IPCV], dx1[IPCV], dx2[IPCV];
	double ay[IPCV - 1], by[IPCV - 1], cy[IPCV - 1], dy[IPCV - 1], d1y[IPCV - 1], d2y[IPCV - 1];
	double axy[JPCV], bxy[JPCV], cxy[JPCV], dxy[JPCV], d1xy[JPCV], d2xy[JPCV];
	int itn = 0;
	int itall = 0;

	////////Regulate of diffusion term of u////////////
	for (i = 0; i < IPCV; i++)
	{
		for (j = 0; j < JPCV - 1; j++)
		{
			dw[i][j] = dely / delx;
			de[i][j] = dely / delx;
			ds[i][j] = delx / dely;
			dn[i][j] = delx / dely;
		}
	}
	i = 0;
	for (j = 1; j <= JPCV - 3; j++)
	{
		ds[i][j] = 2 * delx / dely;
	}
	///************////////¡ü¡üset the 1*1 CV on bottom
	i = IPCV - 1;
	for (j = 1; j <= JPCV - 3; j++)
	{
		dn[i][j] = 2 * delx / dely;
	}
	///************//////// ¡ü¡üset the 1*1 CV on top
	j = 0;
	for (i = 0; i <= IPCV - 1; i++)
	{
		ds[i][j] = 1.5*delx / dely;
		dn[i][j] = 1.5*delx / dely;
	}
	///************//////// ¡ü¡üset the 1*1.5 CV on left
	j = JPCV - 2;
	for (i = 0; i <= IPCV - 1; i++)
	{
		ds[i][j] = 1.5*delx / dely;
		dn[i][j] = 1.5*delx / dely;
	}
	///************//////// ¡ü¡üset the 1*1.5 CV on right
	ds[0][0] = 3 * delx / dely, dn[IPCV - 1][0] = 3 * delx / dely;
	dn[IPCV - 1][JPCV - 2] = 3 * delx / dely, ds[0][JPCV - 2] = 3 * delx / dely;
	///************//////// ¡ü¡üset the four 1*1.5 CV on corners



	/////////////v efficienct///////////////////////////////////
	double fvw[IPCV - 1][JPCV], fve[IPCV - 1][JPCV], fvs[IPCV - 1][JPCV], fvn[IPCV - 1][JPCV];
	double Pvs[IPCV - 1][JPCV], Pvn[IPCV - 1][JPCV], Pvw[IPCV - 1][JPCV], Pve[IPCV - 1][JPCV];
	double APvs[IPCV - 1][JPCV], APvn[IPCV - 1][JPCV], APvw[IPCV - 1][JPCV], APve[IPCV - 1][JPCV];
	double dvw[IPCV - 1][JPCV], dve[IPCV - 1][JPCV], dvs[IPCV - 1][JPCV], dvn[IPCV - 1][JPCV];
	////////Regulate of diffusion term of v////////////
	for (j = 0; j < JPCV; j++)
	{
		for (i = 0; i < IPCV - 1; i++)
		{
			dvw[i][j] = dely / delx;
			dve[i][j] = dely / delx;
			dvs[i][j] = delx / dely;
			dvn[i][j] = delx / dely;
		}
	}
	j = 0;
	for (i = 1; i <= IPCV - 3; i++)
	{
		dvw[i][j] = 2 * dely / delx;
	}
	///************////////¡ü¡üset the 1*1 CV on left
	j = JPCV - 1;
	for (i = 1; i <= IPCV - 3; i++)
	{
		dve[i][j] = 2 * dely / delx;
	}
	///************//////// ¡ü¡üset the 1*1 CV on right
	i = 0;
	for (j = 0; j <= JPCV - 1; j++)
	{
		dvw[i][j] = 1.5*dely / delx;
		dve[i][j] = 1.5*dely / delx;
	}
	///************//////// ¡ü¡üset the 1*1.5 CV on bottom
	i = IPCV - 2;
	for (j = 0; j <= JPCV - 1; j++)
	{
		dvw[i][j] = 1.5*dely / delx;
		dve[i][j] = 1.5*dely / delx;
	}
	///************//////// ¡ü¡üset the 1*1.5 CV on right
	dvw[0][0] = 3 * dely / delx, dvw[IPCV - 2][0] = 3 * dely / delx;
	dve[IPCV - 2][JPCV - 1] = 3 * dely / delx;
	dve[0][JPCV - 1] = 3 * dely / delx;
	///************//////// ¡ü¡üset the four 1*1.5 CV on corners
	//n-1*N MATRIX


	/////////start calculate T/////////////////
	double dtw[IPCV][JPCV] = { 0 }, dte[IPCV][JPCV] = { 0 }, dts[IPCV][JPCV] = { 0 }, dtn[IPCV][JPCV] = { 0 };
	//regulate these coefficient/////////////
	for (i = 0; i < IPCV; i++)
	{
		for (j = 0; j < JPCV; j++)
		{
			dtw[i][j] = dely / delx;
			dte[i][j] = dely / delx;
			dts[i][j] = delx / dely;
			dtn[i][j] = delx / dely;
		}

	}
	j = 0;
	for (i = 0; i < IPCV; i++)
	{
		dtw[i][j] = 2 * dely / delx;
	}
	j = JPCV - 1;
	for (i = 0; i < IPCV; i++)
	{
		dte[i][j] = 2 * dely / delx;
	}
	i = 0;
	for (j = 0; j < JPCV; j++)
	{
		dts[i][j] = 2 * delx / dely;
	}
	i = IPCV - 1;
	for (j = 0; j < JPCV; j++)
	{
		dtn[i][j] = 2 * delx / dely;
	}
	//////////////////////regulate over////////////////


	ofstream in; in.open("u4.csv", ios::ate);
	////////**********************/the over all iteration start here//////////////
	////////**********************/the over all iteration start here//////////////
	////////**********************/the over all iteration start here//////////////
	///////this program will sweep u and v in both directions, but p only in one direction/////////////
	do
	{
		itall = itall + 1;
		Rp = 0;
		Ru = 0;
		Rv = 0;
		///////////******************start solving momentum of u*******************//////////////
		///////////******************start solving momentum of u*******************//////////////
		///////////******************start solving momentum of u*******************//////////////
		//******************from bottom to top**************//////
		//***********line by line solution**************/////////
		//for (i=1;i<IPCV+1;i++)				//bottom to top
		for (i = IPCV; i > 0; i--)                    //top to bottom
		{
			//left 1.5 CV coefficient//
			j = 1;
			fw[i - 1][j - 1] = rho * u[i][j - 1] * dely;
			Pw[i - 1][j - 1] = fw[i - 1][j - 1] / (mu*dw[i - 1][j - 1]);
			APw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pw[i - 1][j - 1]), 5));
			a[j - 1] = (mu*dw[i - 1][j - 1])*APw[i - 1][j - 1] + max(fw[i - 1][j - 1], 0);
			a[j - 1] = -a[j - 1];

			fe[i - 1][j - 1] = rho * (u[i][j] + u[i][j + 1]) / 2 * dely;
			Pe[i - 1][j - 1] = fe[i - 1][j - 1] / (mu*de[i - 1][j - 1]);
			APe[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pe[i - 1][j - 1]), 5));
			c[j - 1] = (mu*de[i - 1][j - 1])*APe[i - 1][j - 1] + max(-fe[i - 1][j - 1], 0);
			c[j - 1] = -c[j - 1];

			fs[i - 1][j - 1] = rho * (v[i - 1][j] + v[i - 1][j + 1]) / 2 * delx + rho * (v[i - 1][j - 1] + v[i - 1][j]) / 2 * delx / 2;
			Ps[i - 1][j - 1] = fs[i - 1][j - 1] / (mu*ds[i - 1][j - 1]);
			APs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ps[i - 1][j - 1]), 5));
			d1[j - 1] = (mu*ds[i - 1][j - 1])*APs[i - 1][j - 1] + max(fs[i - 1][j - 1], 0);

			fn[i - 1][j - 1] = rho * (v[i][j] + v[i][j + 1]) / 2 * delx + rho * (v[i][j - 1] + v[i][j]) / 2 * delx / 2;
			Pn[i - 1][j - 1] = fn[i - 1][j - 1] / (mu*dn[i - 1][j - 1]);
			APn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pn[i - 1][j - 1]), 5));
			d2[j - 1] = (mu*dn[i - 1][j - 1])*APn[i - 1][j - 1] + max(-fn[i - 1][j - 1], 0);

			b[j - 1] = (-a[j - 1] - c[j - 1] + d1[j - 1] + d2[j - 1]) / w;        //underrelax
			d[j - 1] = d1[j - 1] * u[i - 1][j] + d2[j - 1] * u[i + 1][j] + (-a[j - 1])*u[i][j - 1] + (p[i][j] - p[i][j + 1])*dely + b[j - 1] * w*u[i][j] * (1 / w - 1);


			//////middle 1 CV coefficient//////////

			for (j = 2; j < JPCV - 1; j++)
			{
				fw[i - 1][j - 1] = rho * (u[i][j] + u[i][j - 1]) / 2 * dely;
				Pw[i - 1][j - 1] = fw[i - 1][j - 1] / (mu*dw[i - 1][j - 1]);
				APw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pw[i - 1][j - 1]), 5));
				a[j - 1] = (mu*dw[i - 1][j - 1])*APw[i - 1][j - 1] + max(fw[i - 1][j - 1], 0);
				a[j - 1] = -a[j - 1];

				fe[i - 1][j - 1] = rho * (u[i][j] + u[i][j + 1]) / 2 * dely;
				Pe[i - 1][j - 1] = fe[i - 1][j - 1] / (mu*de[i - 1][j - 1]);
				APe[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pe[i - 1][j - 1]), 5));
				c[j - 1] = (mu*de[i - 1][j - 1])*APe[i - 1][j - 1] + max(-fe[i - 1][j - 1], 0);
				c[j - 1] = -c[j - 1];

				fs[i - 1][j - 1] = rho * (v[i - 1][j] + v[i - 1][j + 1]) / 2 * delx;
				Ps[i - 1][j - 1] = fs[i - 1][j - 1] / (mu*ds[i - 1][j - 1]);
				APs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ps[i - 1][j - 1]), 5));
				d1[j - 1] = (mu*ds[i - 1][j - 1])*APs[i - 1][j - 1] + max(fs[i - 1][j - 1], 0);

				fn[i - 1][j - 1] = rho * (v[i][j] + v[i][j + 1]) / 2 * delx;
				Pn[i - 1][j - 1] = fn[i - 1][j - 1] / (mu*dn[i - 1][j - 1]);
				APn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pn[i - 1][j - 1]), 5));
				d2[j - 1] = (mu*dn[i - 1][j - 1])*APn[i - 1][j - 1] + max(-fn[i - 1][j - 1], 0);

				b[j - 1] = (-a[j - 1] - c[j - 1] + d1[j - 1] + d2[j - 1]) / w;        //underrelax
				d[j - 1] = d1[j - 1] * u[i - 1][j] + d2[j - 1] * u[i + 1][j] + (p[i][j] - p[i][j + 1])*dely + b[j - 1] * w*u[i][j] * (1 / w - 1);
			}
			////middle ones end///

			//right 1.5CV Coefficient//
			j = JPCV - 1;
			fw[i - 1][j - 1] = rho * (u[i][j] + u[i][j - 1]) / 2 * dely;
			Pw[i - 1][j - 1] = fw[i - 1][j - 1] / (mu*dw[i - 1][j - 1]);
			APw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pw[i - 1][j - 1]), 5));
			a[j - 1] = (mu*dw[i - 1][j - 1])*APw[i - 1][j - 1] + max(fw[i - 1][j - 1], 0);
			a[j - 1] = -a[j - 1];

			fe[i - 1][j - 1] = rho * u[i][j + 1] * dely;
			Pe[i - 1][j - 1] = fe[i - 1][j - 1] / (mu*de[i - 1][j - 1]);
			APe[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pe[i - 1][j - 1]), 5));
			c[j - 1] = (mu*de[i - 1][j - 1])*APe[i - 1][j - 1] + max(-fe[i - 1][j - 1], 0);
			c[j - 1] = 0;       //outflow


			fs[i - 1][j - 1] = rho * (v[i - 1][j] + v[i - 1][j + 1]) / 2 * delx + rho * (v[i - 1][j + 1] + v[i - 1][j + 2]) / 2 / 2 * delx;
			Ps[i - 1][j - 1] = fs[i - 1][j - 1] / (mu*ds[i - 1][j - 1]);
			APs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ps[i - 1][j - 1]), 5));
			d1[j - 1] = (mu*ds[i - 1][j - 1])*APs[i - 1][j - 1] + max(fs[i - 1][j - 1], 0);

			fn[i - 1][j - 1] = rho * (v[i][j] + v[i][j + 1]) / 2 * delx + rho * (v[i][j + 1] + v[i][j + 2]) / 2 / 2 * delx;
			Pn[i - 1][j - 1] = fn[i - 1][j - 1] / (mu*dn[i - 1][j - 1]);
			APn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pn[i - 1][j - 1]), 5));
			d2[j - 1] = (mu*dn[i - 1][j - 1])*APn[i - 1][j - 1] + max(-fn[i - 1][j - 1], 0);

			b[j - 1] = (-a[j - 1] - c[j - 1] + d1[j - 1] + d2[j - 1]) / w;        //underrelax
			d[j - 1] = d1[j - 1] * u[i - 1][j] + d2[j - 1] * u[i + 1][j] + (-c[j - 1])*u[i][j + 1] + (p[i][j] - p[i][j + 1])*dely + b[j - 1] * w*u[i][j] * (1 / w - 1);

			////////////////////set du ready/////////////////
			for (id = 1; id < JPCV; id++)
			{
				du[i][id] = dely / (b[id - 1] * w);
			}
			////////////////////set du ready/////////////////

			//////***********coefficient end******begin solve///////	
			double *x = new double[JPCV - 1];
			TDMAJ1(x, a, b, c, d);
			for (k = 0; k < JPCV - 1; k++)
			{
				u[i][k + 1] = x[k];
			}
			delete[] x;
		}////////**botton to top end********///





		//////////////Left to right/////////////
		////-as*Ts+ap/w*Tp-an*Tn=aw*Tw+ae*Te+b+ap*Tp*(1/w-1)  
		j = 1;
		for (i = 1; i < IPCV + 1; i++)
		{
			fw[i - 1][j - 1] = rho * u[i][j - 1] * dely;
			Pw[i - 1][j - 1] = fw[i - 1][j - 1] / (mu*dw[i - 1][j - 1]);
			APw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pw[i - 1][j - 1]), 5));
			dx1[i - 1] = (mu*dw[i - 1][j - 1])*APw[i - 1][j - 1] + max(fw[i - 1][j - 1], 0);

			fe[i - 1][j - 1] = rho * (u[i][j] + u[i][j + 1]) / 2 * dely;
			Pe[i - 1][j - 1] = fe[i - 1][j - 1] / (mu*de[i - 1][j - 1]);
			APe[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pe[i - 1][j - 1]), 5));
			dx2[i - 1] = (mu*de[i - 1][j - 1])*APe[i - 1][j - 1] + max(-fe[i - 1][j - 1], 0);

			fs[i - 1][j - 1] = rho * (v[i - 1][j] + v[i - 1][j + 1]) / 2 * delx + rho * (v[i - 1][j - 1] + v[i - 1][j]) / 2 * delx / 2;
			Ps[i - 1][j - 1] = fs[i - 1][j - 1] / (mu*ds[i - 1][j - 1]);
			APs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ps[i - 1][j - 1]), 5));
			ax[i - 1] = (mu*ds[i - 1][j - 1])*APs[i - 1][j - 1] + max(fs[i - 1][j - 1], 0);
			ax[i - 1] = -ax[i - 1];

			fn[i - 1][j - 1] = rho * (v[i][j] + v[i][j + 1]) / 2 * delx + rho * (v[i][j - 1] + v[i][j]) / 2 * delx / 2;
			Pn[i - 1][j - 1] = fn[i - 1][j - 1] / (mu*dn[i - 1][j - 1]);
			APn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pn[i - 1][j - 1]), 5));
			cx[i - 1] = (mu*dn[i - 1][j - 1])*APn[i - 1][j - 1] + max(-fn[i - 1][j - 1], 0);
			cx[i - 1] = -cx[i - 1];

			bx[i - 1] = (-ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1]) / w;        //underrelax
			dx[i - 1] = dx1[i - 1] * u[i][j - 1] + dx2[i - 1] * u[i][j + 1] + (p[i][j] - p[i][j + 1])*dely + bx[i - 1] * w*u[i][j] * (1 / w - 1);
		}
		i = 1;
		dx[i - 1] = dx1[i - 1] * u[i][j - 1] + dx2[i - 1] * u[i][j + 1] + (-ax[i - 1])*u[i - 1][j] + (p[i][j] - p[i][j + 1])*dely + bx[i - 1] * w*u[i][j] * (1 / w - 1);
		i = IPCV;
		dx2[i - 1] = 0;
		bx[i - 1] = (-ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1]) / w;        //underrelax
		dx[i - 1] = dx1[i - 1] * u[i][j - 1] + dx2[i - 1] * u[i][j + 1] + (-cx[i - 1])*u[i + 1][j] + (p[i][j] - p[i][j + 1])*dely + bx[i - 1] * w*u[i][j] * (1 / w - 1);
		//////***********coefficient end******begin solve///////
		double *y = new double[IPCV];
		TDMAI(y, ax, bx, cx, dx);
		for (k = 0; k < IPCV; k++)
		{
			u[k + 1][j] = y[k];
		}
		///////***********************left over*********************/////////////////////
		for (j = 2; j < JPCV - 1; j++)
		{
			for (i = 1; i < IPCV + 1; i++)
			{
				fw[i - 1][j - 1] = rho * (u[i][j] + u[i][j - 1]) / 2 * dely;
				Pw[i - 1][j - 1] = fw[i - 1][j - 1] / (mu*dw[i - 1][j - 1]);
				APw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pw[i - 1][j - 1]), 5));
				dx1[i - 1] = (mu*dw[i - 1][j - 1])*APw[i - 1][j - 1] + max(fw[i - 1][j - 1], 0);


				fe[i - 1][j - 1] = rho * (u[i][j] + u[i][j + 1]) / 2 * dely;
				Pe[i - 1][j - 1] = fe[i - 1][j - 1] / (mu*de[i - 1][j - 1]);
				APe[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pe[i - 1][j - 1]), 5));
				dx2[i - 1] = (mu*de[i - 1][j - 1])*APe[i - 1][j - 1] + max(-fe[i - 1][j - 1], 0);


				fs[i - 1][j - 1] = rho * (v[i - 1][j] + v[i - 1][j + 1]) / 2 * delx;
				Ps[i - 1][j - 1] = fs[i - 1][j - 1] / (mu*ds[i - 1][j - 1]);
				APs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ps[i - 1][j - 1]), 5));
				ax[i - 1] = (mu*ds[i - 1][j - 1])*APs[i - 1][j - 1] + max(fs[i - 1][j - 1], 0);
				ax[i - 1] = -ax[i - 1];

				fn[i - 1][j - 1] = rho * (v[i][j] + v[i][j + 1]) / 2 * delx;
				Pn[i - 1][j - 1] = fn[i - 1][j - 1] / (mu*dn[i - 1][j - 1]);
				APn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pn[i - 1][j - 1]), 5));
				cx[i - 1] = (mu*dn[i - 1][j - 1])*APn[i - 1][j - 1] + max(-fn[i - 1][j - 1], 0);
				cx[i - 1] = -cx[i - 1];

				bx[i - 1] = (-ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1]) / w;        //underrelax
				dx[i - 1] = dx1[i - 1] * u[i][j - 1] + dx2[i - 1] * u[i][j + 1] + (p[i][j] - p[i][j + 1])*dely + bx[i - 1] * w*u[i][j] * (1 / w - 1);
			}
			i = 1;
			dx[i - 1] = dx1[i - 1] * u[i][j - 1] + dx2[i - 1] * u[i][j + 1] + (-ax[i - 1])*u[i - 1][j] + (p[i][j] - p[i][j + 1])*dely + bx[i - 1] * w*u[i][j] * (1 / w - 1);
			i = IPCV;
			dx2[i - 1] = 0;
			bx[i - 1] = (-ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1]) / w;        //underrelax
			dx[i - 1] = dx1[i - 1] * u[i][j - 1] + dx2[i - 1] * u[i][j + 1] + (-cx[i - 1])*u[i + 1][j] + (p[i][j] - p[i][j + 1])*dely + bx[i - 1] * w*u[i][j] * (1 / w - 1);
			//////***********coefficient end******begin solve///////
			TDMAI(y, ax, bx, cx, dx);
			for (k = 0; k < IPCV; k++)
			{
				u[k + 1][j] = y[k];
			}
		}

		/////**************middle over**************/////////////////////////////////
		j = JPCV - 1;
		for (i = 1; i < IPCV + 1; i++)
		{
			fw[i - 1][j - 1] = rho * (u[i][j] + u[i][j - 1]) / 2 * dely;
			Pw[i - 1][j - 1] = fw[i - 1][j - 1] / (mu*dw[i - 1][j - 1]);
			APw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pw[i - 1][j - 1]), 5));
			dx1[i - 1] = (mu*dw[i - 1][j - 1])*APw[i - 1][j - 1] + max(fw[i - 1][j - 1], 0);

			fe[i - 1][j - 1] = rho * u[i][j + 1] * dely;
			Pe[i - 1][j - 1] = fe[i - 1][j - 1] / (mu*de[i - 1][j - 1]);
			APe[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pe[i - 1][j - 1]), 5));
			dx2[i - 1] = (mu*de[i - 1][j - 1])*APe[i - 1][j - 1] + max(-fe[i - 1][j - 1], 0);

			fs[i - 1][j - 1] = rho * (v[i - 1][j] + v[i - 1][j + 1]) / 2 * delx + rho * (v[i - 1][j + 1] + v[i - 1][j + 2]) / 2 / 2 * delx;
			Ps[i - 1][j - 1] = fs[i - 1][j - 1] / (mu*ds[i - 1][j - 1]);
			APs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ps[i - 1][j - 1]), 5));
			ax[i - 1] = (mu*ds[i - 1][j - 1])*APs[i - 1][j - 1] + max(fs[i - 1][j - 1], 0);
			ax[i - 1] = -ax[i - 1];

			fn[i - 1][j - 1] = rho * (v[i][j] + v[i][j + 1]) / 2 * delx + rho * (v[i][j + 1] + v[i][j + 2]) / 2 / 2 * delx;
			Pn[i - 1][j - 1] = fn[i - 1][j - 1] / (mu*dn[i - 1][j - 1]);
			APn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pn[i - 1][j - 1]), 5));
			cx[i - 1] = (mu*dn[i - 1][j - 1])*APn[i - 1][j - 1] + max(-fn[i - 1][j - 1], 0);
			cx[i - 1] = -cx[i - 1];

			bx[i - 1] = (-ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1]) / w;        //underrelax
			dx[i - 1] = dx1[i - 1] * u[i][j - 1] + dx2[i - 1] * u[i][j + 1] + (p[i][j] - p[i][j + 1])*dely + bx[i - 1] * w*u[i][j] * (1 / w - 1);
		}
		i = 1;
		dx[i - 1] = dx1[i - 1] * u[i][j - 1] + dx2[i - 1] * u[i][j + 1] + (-ax[i - 1])*u[i - 1][j] + (p[i][j] - p[i][j + 1])*dely + bx[i - 1] * w*u[i][j] * (1 / w - 1);
		i = IPCV;
		dx2[i - 1] = 0;
		bx[i - 1] = (-ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1]) / w;        //underrelax
		dx[i - 1] = dx1[i - 1] * u[i][j - 1] + dx2[i - 1] * u[i][j + 1] + (-cx[i - 1])*u[i + 1][j] + (p[i][j] - p[i][j + 1])*dely + bx[i - 1] * w*u[i][j] * (1 / w - 1);
		//////***********coefficient end******begin solve///////
		TDMAI(y, ax, bx, cx, dx);
		for (k = 0; k < IPCV; k++)
		{
			u[k + 1][j] = y[k];
		}
		//////////////////////right side over/////////////////////////////////////
		///////////////outflow correction///////////////////////////////////////
		double massratio, flowin = 0, flowout = 0;
		flowin = IPCV * Vinf;
		for (i = 0; i < IPCV + 2; i++)
		{
			flowout = flowout + u[i][JPCV - 1];
		}
		massratio = flowin / (flowout + 1e-15);
		for (i = 0; i < IPCV + 2; i++)
		{
			u[i][JPCV] = u[i][JPCV - 1] * massratio;
		}
		///////////////outflow correction///////////////////////////////////////					

		/////////////////calculate th convergence of u////////////////////
		for (i = 1; i < IPCV + 1; i++)
		{
			for (j = 1; j < JPCV + 1; j++)
			{
				uold[i][j] = u[i][j];
			}
		}


		///////////////////**************solving v momentum/***********************************///////////////
		///////////////////**************solving v momentum/***********************************///////////////
		///////////////////**************solving v momentum/***********************************///////////////
		///////////////////**************solving v momentum/***********************************///////////////
		///////////////////**************solving v momentum/***********************************///////////////
		///////////////////**************solving v momentum/***********************************///////////////

		///////////////////sweep from left to right//////////////////////////////////
		///////////////////I will sweep twice in this direction//////////////////////
		//-as*Ts+ap/w*Tp-an*Tn=aw*Tw+ae*Te+b+ap*Tp*(1/w-1)    

		for (j = 1; j < JPCV; j++)
		{
			///////////regulate the bottom//////////////////
			i = 1;
			fvw[i - 1][j - 1] = rho * (u[i + 1][j - 1] + u[i][j - 1]) / 2 * dely + rho * (u[i][j - 1] + u[i - 1][j - 1]) / 2 / 2 * dely;
			Pvw[i - 1][j - 1] = fvw[i - 1][j - 1] / (mu*dvw[i - 1][j - 1]);
			APvw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvw[i - 1][j - 1]), 5));
			d1y[i - 1] = (mu*dvw[i - 1][j - 1])*APvw[i - 1][j - 1] + max(fvw[i - 1][j - 1], 0);


			fve[i - 1][j - 1] = rho * (u[i + 1][j] + u[i][j]) / 2 * dely + rho * (u[i][j] + u[i - 1][j]) / 2 / 2 * dely;
			Pve[i - 1][j - 1] = fve[i - 1][j - 1] / (mu*dve[i - 1][j - 1]);
			APve[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pve[i - 1][j - 1]), 5));
			d2y[i - 1] = (mu*dve[i - 1][j - 1])*APve[i - 1][j - 1] + max(-fve[i - 1][j - 1], 0);


			fvs[i - 1][j - 1] = rho * v[i - 1][j] * delx;
			Pvs[i - 1][j - 1] = fvs[i - 1][j - 1] / (mu*dvs[i - 1][j - 1]);
			APvs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvs[i - 1][j - 1]), 5));
			ay[i - 1] = (mu*dvs[i - 1][j - 1])*APvs[i - 1][j - 1] + max(fvs[i - 1][j - 1], 0);
			ay[i - 1] = -ay[i - 1];

			fvn[i - 1][j - 1] = rho * (v[i][j] + v[i + 1][j]) / 2 * delx;
			Pvn[i - 1][j - 1] = fvn[i - 1][j - 1] / (mu*dvn[i - 1][j - 1]);
			APvn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvn[i - 1][j - 1]), 5));
			cy[i - 1] = (mu*dvn[i - 1][j - 1])*APvn[i - 1][j - 1] + max(-fvn[i - 1][j - 1], 0);
			cy[i - 1] = -cy[i - 1];

			by[i - 1] = (-ay[i - 1] - cy[i - 1] + d1y[i - 1] + d2y[i - 1]) / w;
			dy[i - 1] = d1y[i - 1] * v[i][j - 1] + d2y[i - 1] * v[i][j + 1] + (-ay[i - 1])*v[i - 1][j] + (p[i][j] - p[i + 1][j])*delx + by[i - 1] * w*v[i][j] * (1 / w - 1);





			///////////regulate the middle//////////////////
			for (i = 2; i < IPCV - 1; i++)
			{
				fvw[i - 1][j - 1] = rho * (u[i][j - 1] + u[i + 1][j - 1]) / 2 * dely;
				Pvw[i - 1][j - 1] = fvw[i - 1][j - 1] / (mu*dvw[i - 1][j - 1]);
				APvw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvw[i - 1][j - 1]), 5));
				d1y[i - 1] = (mu*dvw[i - 1][j - 1])*APvw[i - 1][j - 1] + max(fvw[i - 1][j - 1], 0);

				fve[i - 1][j - 1] = rho * (u[i][j] + u[i + 1][j]) / 2 * dely;
				Pve[i - 1][j - 1] = fve[i - 1][j - 1] / (mu*dve[i - 1][j - 1]);
				APve[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pve[i - 1][j - 1]), 5));
				d2y[i - 1] = (mu*dve[i - 1][j - 1])*APve[i - 1][j - 1] + max(-fve[i - 1][j - 1], 0);

				fvs[i - 1][j - 1] = rho * (v[i][j] + v[i - 1][j]) / 2 * delx;
				Pvs[i - 1][j - 1] = fvs[i - 1][j - 1] / (mu*dvs[i - 1][j - 1]);
				APvs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvs[i - 1][j - 1]), 5));
				ay[i - 1] = (mu*dvs[i - 1][j - 1])*APvs[i - 1][j - 1] + max(fvs[i - 1][j - 1], 0);
				ay[i - 1] = -ay[i - 1];

				fvn[i - 1][j - 1] = rho * (v[i][j] + v[i + 1][j]) / 2 * delx;
				Pvn[i - 1][j - 1] = fvn[i - 1][j - 1] / (mu*dvn[i - 1][j - 1]);
				APvn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvn[i - 1][j - 1]), 5));
				cy[i - 1] = (mu*dvn[i - 1][j - 1])*APvn[i - 1][j - 1] + max(-fvn[i - 1][j - 1], 0);
				cy[i - 1] = -cy[i - 1];

				by[i - 1] = (-ay[i - 1] - cy[i - 1] + d1y[i - 1] + d2y[i - 1]) / w;
				dy[i - 1] = d1y[i - 1] * v[i][j - 1] + d2y[i - 1] * v[i][j + 1] + (p[i][j] - p[i + 1][j])*delx + by[i - 1] * w*v[i][j] * (1 / w - 1);


			}
			///////////regulate the top//////////////////
			i = IPCV - 1;
			fvw[i - 1][j - 1] = rho * (u[i][j - 1] + u[i + 1][j - 1]) / 2 * dely + rho * (u[i + 1][j - 1] + u[i + 2][j - 1]) / 2 / 2 * dely;
			Pvw[i - 1][j - 1] = fvw[i - 1][j - 1] / (mu*dvw[i - 1][j - 1]);
			APvw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvw[i - 1][j - 1]), 5));
			d1y[i - 1] = (mu*dvw[i - 1][j - 1])*APvw[i - 1][j - 1] + max(fvw[i - 1][j - 1], 0);


			fve[i - 1][j - 1] = rho * (u[i][j] + u[i + 1][j]) / 2 * dely + rho * (u[i + 1][j] + u[i + 2][j]) / 2 / 2 * dely;
			Pve[i - 1][j - 1] = fve[i - 1][j - 1] / (mu*dve[i - 1][j - 1]);
			APve[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pve[i - 1][j - 1]), 5));
			d2y[i - 1] = (mu*dve[i - 1][j - 1])*APve[i - 1][j - 1] + max(-fve[i - 1][j - 1], 0);


			fvs[i - 1][j - 1] = rho * (v[i][j] + v[i - 1][j]) / 2 * delx;
			Pvs[i - 1][j - 1] = fvs[i - 1][j - 1] / (mu*dvs[i - 1][j - 1]);
			APvs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvs[i - 1][j - 1]), 5));
			ay[i - 1] = (mu*dvs[i - 1][j - 1])*APvs[i - 1][j - 1] + max(fvs[i - 1][j - 1], 0);
			ay[i - 1] = -ay[i - 1];

			fvn[i - 1][j - 1] = rho * v[i + 1][j] * delx;
			Pvn[i - 1][j - 1] = fvn[i - 1][j - 1] / (mu*dvn[i - 1][j - 1]);
			APvn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvn[i - 1][j - 1]), 5));
			cy[i - 1] = (mu*dvn[i - 1][j - 1])*APvn[i - 1][j - 1] + max(-fvn[i - 1][j - 1], 0);
			cy[i - 1] = -cy[i - 1];

			by[i - 1] = (-ay[i - 1] - cy[i - 1] + d1y[i - 1] + d2y[i - 1]) / w;
			dy[i - 1] = d1y[i - 1] * v[i][j - 1] + d2y[i - 1] * v[i][j + 1] + (-cy[i - 1])*v[i + 1][j] + (p[i][j] - p[i + 1][j])*delx + by[i - 1] * w*v[i][j] * (1 / w - 1);

			/////////////****set dv ready////////////
			for (id = 1; id < IPCV; id++)
			{
				dv[id][j] = delx / (by[id - 1] * w);
			}
			/////////////****set dv ready////////////

			//////////////*********solve it***************//////////////
			double *z = new double[IPCV - 1];
			TDMAI1(z, ay, by, cy, dy);            //****************need change
			for (k = 0; k < IPCV - 1; k++)
			{
				v[k + 1][j] = z[k];
			}
			delete[] z;

		}

		j = JPCV;
		///////////regulate the bottom//////////////////
		i = 1;
		fvw[i - 1][j - 1] = rho * (u[i + 1][j - 1] + u[i][j - 1]) / 2 * dely + rho * (u[i][j - 1] + u[i - 1][j - 1]) / 2 / 2 * dely;
		Pvw[i - 1][j - 1] = fvw[i - 1][j - 1] / (mu*dvw[i - 1][j - 1]);
		APvw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvw[i - 1][j - 1]), 5));
		d1y[i - 1] = (mu*dvw[i - 1][j - 1])*APvw[i - 1][j - 1] + max(fvw[i - 1][j - 1], 0);


		fve[i - 1][j - 1] = rho * (u[i + 1][j] + u[i][j]) / 2 * dely + rho * (u[i][j] + u[i - 1][j]) / 2 / 2 * dely;
		Pve[i - 1][j - 1] = fve[i - 1][j - 1] / (mu*dve[i - 1][j - 1]);
		APve[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pve[i - 1][j - 1]), 5));
		d2y[i - 1] = 0;


		fvs[i - 1][j - 1] = rho * v[i - 1][j] * delx;
		Pvs[i - 1][j - 1] = fvs[i - 1][j - 1] / (mu*dvs[i - 1][j - 1]);
		APvs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvs[i - 1][j - 1]), 5));
		ay[i - 1] = (mu*dvs[i - 1][j - 1])*APvs[i - 1][j - 1] + max(fvs[i - 1][j - 1], 0);
		ay[i - 1] = -ay[i - 1];

		fvn[i - 1][j - 1] = rho * (v[i][j] + v[i + 1][j]) / 2 * delx;
		Pvn[i - 1][j - 1] = fvn[i - 1][j - 1] / (mu*dvn[i - 1][j - 1]);
		APvn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvn[i - 1][j - 1]), 5));
		cy[i - 1] = (mu*dvn[i - 1][j - 1])*APvn[i - 1][j - 1] + max(-fvn[i - 1][j - 1], 0);
		cy[i - 1] = -cy[i - 1];

		by[i - 1] = (-ay[i - 1] - cy[i - 1] + d1y[i - 1] + d2y[i - 1]) / w;
		dy[i - 1] = d1y[i - 1] * v[i][j - 1] + d2y[i - 1] * v[i][j + 1] + (-ay[i - 1])*v[i - 1][j] + (p[i][j] - p[i + 1][j])*delx + by[i - 1] * w*v[i][j] * (1 / w - 1);
		///////////regulate the middle//////////////////
		for (i = 2; i < IPCV - 1; i++)
		{
			fvw[i - 1][j - 1] = rho * (u[i][j - 1] + u[i + 1][j - 1]) / 2 * dely;
			Pvw[i - 1][j - 1] = fvw[i - 1][j - 1] / (mu*dvw[i - 1][j - 1]);
			APvw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvw[i - 1][j - 1]), 5));
			d1y[i - 1] = (mu*dvw[i - 1][j - 1])*APvw[i - 1][j - 1] + max(fvw[i - 1][j - 1], 0);

			fve[i - 1][j - 1] = rho * (u[i][j] + u[i + 1][j]) / 2 * dely;
			Pve[i - 1][j - 1] = fve[i - 1][j - 1] / (mu*dve[i - 1][j - 1]);
			APve[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pve[i - 1][j - 1]), 5));
			d2y[i - 1] = 0;

			fvs[i - 1][j - 1] = rho * (v[i][j] + v[i - 1][j]) / 2 * delx;
			Pvs[i - 1][j - 1] = fvs[i - 1][j - 1] / (mu*dvs[i - 1][j - 1]);
			APvs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvs[i - 1][j - 1]), 5));
			ay[i - 1] = (mu*dvs[i - 1][j - 1])*APvs[i - 1][j - 1] + max(fvs[i - 1][j - 1], 0);
			ay[i - 1] = -ay[i - 1];

			fvn[i - 1][j - 1] = rho * (v[i][j] + v[i + 1][j]) / 2 * delx;
			Pvn[i - 1][j - 1] = fvn[i - 1][j - 1] / (mu*dvn[i - 1][j - 1]);
			APvn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvn[i - 1][j - 1]), 5));
			cy[i - 1] = (mu*dvn[i - 1][j - 1])*APvn[i - 1][j - 1] + max(-fvn[i - 1][j - 1], 0);
			cy[i - 1] = -cy[i - 1];

			by[i - 1] = (-ay[i - 1] - cy[i - 1] + d1y[i - 1] + d2y[i - 1]) / w;
			dy[i - 1] = d1y[i - 1] * v[i][j - 1] + d2y[i - 1] * v[i][j + 1] + (p[i][j] - p[i + 1][j])*delx + by[i - 1] * w*v[i][j] * (1 / w - 1);


		}
		///////////regulate the top//////////////////
		i = IPCV - 1;
		fvw[i - 1][j - 1] = rho * (u[i][j - 1] + u[i + 1][j - 1]) / 2 * dely + rho * (u[i + 1][j - 1] + u[i + 2][j - 1]) / 2 / 2 * dely;
		Pvw[i - 1][j - 1] = fvw[i - 1][j - 1] / (mu*dvw[i - 1][j - 1]);
		APvw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvw[i - 1][j - 1]), 5));
		d1y[i - 1] = (mu*dvw[i - 1][j - 1])*APvw[i - 1][j - 1] + max(fvw[i - 1][j - 1], 0);


		fve[i - 1][j - 1] = rho * (u[i][j] + u[i + 1][j]) / 2 * dely + rho * (u[i + 1][j] + u[i + 2][j]) / 2 / 2 * dely;
		Pve[i - 1][j - 1] = fve[i - 1][j - 1] / (mu*dve[i - 1][j - 1]);
		APve[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pve[i - 1][j - 1]), 5));
		d2y[i - 1] = 0;


		fvs[i - 1][j - 1] = rho * (v[i][j] + v[i - 1][j]) / 2 * delx;
		Pvs[i - 1][j - 1] = fvs[i - 1][j - 1] / (mu*dvs[i - 1][j - 1]);
		APvs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvs[i - 1][j - 1]), 5));
		ay[i - 1] = (mu*dvs[i - 1][j - 1])*APvs[i - 1][j - 1] + max(fvs[i - 1][j - 1], 0);
		ay[i - 1] = -ay[i - 1];

		fvn[i - 1][j - 1] = rho * v[i + 1][j] * delx;
		Pvn[i - 1][j - 1] = fvn[i - 1][j - 1] / (mu*dvn[i - 1][j - 1]);
		APvn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvn[i - 1][j - 1]), 5));
		cy[i - 1] = (mu*dvn[i - 1][j - 1])*APvn[i - 1][j - 1] + max(-fvn[i - 1][j - 1], 0);
		cy[i - 1] = -cy[i - 1];

		by[i - 1] = (-ay[i - 1] - cy[i - 1] + d1y[i - 1] + d2y[i - 1]) / w;
		dy[i - 1] = d1y[i - 1] * v[i][j - 1] + d2y[i - 1] * v[i][j + 1] + (-cy[i - 1])*v[i + 1][j] + (p[i][j] - p[i + 1][j])*delx + by[i - 1] * w*v[i][j] * (1 / w - 1);

		/////////////****set dv ready////////////
		for (id = 1; id < IPCV; id++)
		{
			dv[id][j] = delx / (by[id - 1] * w);
		}
		/////////////****set dv ready////////////

		//////////////*********solve it***************//////////////
		double *z = new double[IPCV - 1];
		TDMAI1(z, ay, by, cy, dy);            //****************need change
		for (k = 0; k < IPCV - 1; k++)
		{
			v[k + 1][j] = z[k];
		}
		delete[] z;

		///////////////////sweep from bottom to top//////////////////////////////////
		//-aw*Tw+ap/w*Tp-ae*Te=as*Ts+an*Tn+b+ap*Tp*(1/w-1)    [1]
		i = 1;
		for (j = 1; j < JPCV + 1; j++)
		{
			fvw[i - 1][j - 1] = rho * (u[i + 1][j - 1] + u[i][j - 1]) / 2 * dely + rho * (u[i][j - 1] + u[i - 1][j - 1]) / 2 / 2 * dely;
			Pvw[i - 1][j - 1] = fvw[i - 1][j - 1] / (mu*dvw[i - 1][j - 1]);
			APvw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvw[i - 1][j - 1]), 5));
			axy[j - 1] = (mu*dvw[i - 1][j - 1])*APvw[i - 1][j - 1] + max(fvw[i - 1][j - 1], 0);
			axy[j - 1] = -axy[j - 1];

			fve[i - 1][j - 1] = rho * (u[i + 1][j] + u[i][j]) / 2 * dely + rho * (u[i][j] + u[i - 1][j]) / 2 / 2 * dely;
			Pve[i - 1][j - 1] = fve[i - 1][j - 1] / (mu*dve[i - 1][j - 1]);
			APve[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pve[i - 1][j - 1]), 5));
			cxy[j - 1] = (mu*dve[i - 1][j - 1])*APve[i - 1][j - 1] + max(-fve[i - 1][j - 1], 0);
			cxy[j - 1] = -cxy[j - 1];

			fvs[i - 1][j - 1] = rho * v[i - 1][j] * delx;
			Pvs[i - 1][j - 1] = fvs[i - 1][j - 1] / (mu*dvs[i - 1][j - 1]);
			APvs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvs[i - 1][j - 1]), 5));
			d1xy[j - 1] = (mu*dvs[i - 1][j - 1])*APvs[i - 1][j - 1] + max(fvs[i - 1][j - 1], 0);


			fvn[i - 1][j - 1] = rho * (v[i][j] + v[i + 1][j]) / 2 * delx;
			Pvn[i - 1][j - 1] = fvn[i - 1][j - 1] / (mu*dvn[i - 1][j - 1]);
			APvn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvn[i - 1][j - 1]), 5));
			d2xy[j - 1] = (mu*dvn[i - 1][j - 1])*APvn[i - 1][j - 1] + max(-fvn[i - 1][j - 1], 0);


			bxy[j - 1] = (-axy[j - 1] - cxy[j - 1] + d1xy[j - 1] + d2xy[j - 1]) / w;
			dxy[j - 1] = d1xy[j - 1] * v[i - 1][j] + d2xy[j - 1] * v[i + 1][j] + (p[i][j] - p[i + 1][j])*delx + bxy[j - 1] * w*v[i][j] * (1 / w - 1);
		}
		j = 1;
		dxy[j - 1] = d1xy[j - 1] * v[i - 1][j] + d2xy[j - 1] * v[i + 1][j] + (-axy[j - 1])*v[i][j - 1] + (p[i][j] - p[i + 1][j])*delx + bxy[j - 1] * w*v[i][j] * (1 / w - 1);
		j = JPCV;
		cxy[j - 1] = 0;
		bxy[j - 1] = (-axy[j - 1] - cxy[j - 1] + d1xy[j - 1] + d2xy[j - 1]) / w;
		dxy[j - 1] = d1xy[j - 1] * v[i - 1][j] + d2xy[j - 1] * v[i + 1][j] + (-cxy[j - 1])*v[i][j + 1] + (p[i][j] - p[i + 1][j])*delx + bxy[j - 1] * w*v[i][j] * (1 / w - 1);
		/////////**********solve it**********
		double *xy = new double[JPCV];
		TDMAJ(xy, axy, bxy, cxy, dxy);
		for (k = 0; k < JPCV; k++)
		{
			v[i][k + 1] = xy[k];
		}
		/////////////////////////bottom over/////////////////////////////////////////////////
		for (i = 2; i < IPCV - 1; i++)
		{
			for (j = 1; j < JPCV + 1; j++)
			{
				fvw[i - 1][j - 1] = rho * (u[i][j - 1] + u[i + 1][j - 1]) / 2 * dely;
				Pvw[i - 1][j - 1] = fvw[i - 1][j - 1] / (mu*dvw[i - 1][j - 1]);
				APvw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvw[i - 1][j - 1]), 5));
				axy[j - 1] = (mu*dvw[i - 1][j - 1])*APvw[i - 1][j - 1] + max(fvw[i - 1][j - 1], 0);
				axy[j - 1] = -axy[j - 1];

				fve[i - 1][j - 1] = rho * (u[i][j] + u[i + 1][j]) / 2 * dely;
				Pve[i - 1][j - 1] = fve[i - 1][j - 1] / (mu*dve[i - 1][j - 1]);
				APve[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pve[i - 1][j - 1]), 5));
				cxy[j - 1] = (mu*dve[i - 1][j - 1])*APve[i - 1][j - 1] + max(-fve[i - 1][j - 1], 0);
				cxy[j - 1] = -cxy[j - 1];

				fvs[i - 1][j - 1] = rho * (v[i][j] + v[i - 1][j]) / 2 * delx;
				Pvs[i - 1][j - 1] = fvs[i - 1][j - 1] / (mu*dvs[i - 1][j - 1]);
				APvs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvs[i - 1][j - 1]), 5));
				d1xy[j - 1] = (mu*dvs[i - 1][j - 1])*APvs[i - 1][j - 1] + max(fvs[i - 1][j - 1], 0);

				fvn[i - 1][j - 1] = rho * (v[i][j] + v[i + 1][j]) / 2 * delx;
				Pvn[i - 1][j - 1] = fvn[i - 1][j - 1] / (mu*dvn[i - 1][j - 1]);
				APvn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvn[i - 1][j - 1]), 5));
				d2xy[j - 1] = (mu*dvn[i - 1][j - 1])*APvn[i - 1][j - 1] + max(-fvn[i - 1][j - 1], 0);

				bxy[j - 1] = (-axy[j - 1] - cxy[j - 1] + d1xy[j - 1] + d2xy[j - 1]) / w;
				dxy[j - 1] = d1xy[j - 1] * v[i - 1][j] + d2xy[j - 1] * v[i + 1][j] + (p[i][j] - p[i + 1][j])*delx + bxy[j - 1] * w*v[i][j] * (1 / w - 1);
			}
			j = 1;
			dxy[j - 1] = d1xy[j - 1] * v[i - 1][j] + d2xy[j - 1] * v[i + 1][j] + (-axy[j - 1])*v[i][j - 1] + (p[i][j] - p[i + 1][j])*delx + bxy[j - 1] * w*v[i][j] * (1 / w - 1);
			j = JPCV;
			cxy[j - 1] = 0;
			bxy[j - 1] = (-axy[j - 1] - cxy[j - 1] + d1xy[j - 1] + d2xy[j - 1]) / w;
			dxy[j - 1] = d1xy[j - 1] * v[i - 1][j] + d2xy[j - 1] * v[i + 1][j] + (-cxy[j - 1])*v[i][j + 1] + (p[i][j] - p[i + 1][j])*delx + bxy[j - 1] * w*v[i][j] * (1 / w - 1);
			/////////**********solve it**********
			//double *xy = new double[JPCV];
			TDMAJ(xy, axy, bxy, cxy, dxy);
			for (k = 0; k < JPCV; k++)
			{
				v[i][k + 1] = xy[k];
			}

		}
		///////////////************middle over*****************/////////////////////////////
		i = IPCV - 1;
		for (j = 1; j < JPCV + 1; j++)
		{
			fvw[i - 1][j - 1] = rho * (u[i][j - 1] + u[i + 1][j - 1]) / 2 * dely + rho * (u[i + 1][j - 1] + u[i + 2][j - 1]) / 2 / 2 * dely;
			Pvw[i - 1][j - 1] = fvw[i - 1][j - 1] / (mu*dvw[i - 1][j - 1]);
			APvw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvw[i - 1][j - 1]), 5));
			axy[j - 1] = (mu*dvw[i - 1][j - 1])*APvw[i - 1][j - 1] + max(fvw[i - 1][j - 1], 0);
			axy[j - 1] = -axy[j - 1];


			fve[i - 1][j - 1] = rho * (u[i][j] + u[i + 1][j]) / 2 * dely + rho * (u[i + 1][j] + u[i + 2][j]) / 2 / 2 * dely;
			Pve[i - 1][j - 1] = fve[i - 1][j - 1] / (mu*dve[i - 1][j - 1]);
			APve[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pve[i - 1][j - 1]), 5));
			cxy[j - 1] = (mu*dve[i - 1][j - 1])*APve[i - 1][j - 1] + max(-fve[i - 1][j - 1], 0);
			cxy[j - 1] = -cxy[j - 1];


			fvs[i - 1][j - 1] = rho * (v[i][j] + v[i - 1][j]) / 2 * delx;
			Pvs[i - 1][j - 1] = fvs[i - 1][j - 1] / (mu*dvs[i - 1][j - 1]);
			APvs[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvs[i - 1][j - 1]), 5));
			d1xy[j - 1] = (mu*dvs[i - 1][j - 1])*APvs[i - 1][j - 1] + max(fvs[i - 1][j - 1], 0);

			fvn[i - 1][j - 1] = rho * v[i + 1][j] * delx;
			Pvn[i - 1][j - 1] = fvn[i - 1][j - 1] / (mu*dvn[i - 1][j - 1]);
			APvn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pvn[i - 1][j - 1]), 5));
			d2xy[j - 1] = (mu*dvn[i - 1][j - 1])*APvn[i - 1][j - 1] + max(-fvn[i - 1][j - 1], 0);

			bxy[j - 1] = (-axy[j - 1] - cxy[j - 1] + d1xy[j - 1] + d2xy[j - 1]) / w;
			dxy[j - 1] = d1xy[j - 1] * v[i - 1][j] + d2xy[j - 1] * v[i + 1][j] + (p[i][j] - p[i + 1][j])*delx + bxy[j - 1] * w*v[i][j] * (1 / w - 1);
		}
		j = 1;
		dxy[j - 1] = d1xy[j - 1] * v[i - 1][j] + d2xy[j - 1] * v[i + 1][j] + (-axy[j - 1])*v[i][j - 1] + (p[i][j] - p[i + 1][j])*delx + bxy[j - 1] * w*v[i][j] * (1 / w - 1);
		j = JPCV;
		cxy[j - 1] = 0;
		bxy[j - 1] = (-axy[j - 1] - cxy[j - 1] + d1xy[j - 1] + d2xy[j - 1]) / w;
		dxy[j - 1] = d1xy[j - 1] * v[i - 1][j] + d2xy[j - 1] * v[i + 1][j] + (-cxy[j - 1])*v[i][j + 1] + (p[i][j] - p[i + 1][j])*delx + bxy[j - 1] * w*v[i][j] * (1 / w - 1);
		/////////**********solve it**********
		//double *xy = new double[JPCV];
		TDMAJ(xy, axy, bxy, cxy, dxy);
		for (k = 0; k < JPCV; k++)
		{
			v[i][k + 1] = xy[k];
		}
		delete[] xy;
		///////////////////////////top over///////////////////////////////////////////////////////////////////




			/////////////////calculate th convergence of v////////////////////
		for (i = 1; i < IPCV; i++)
		{
			for (j = 1; j < JPCV + 1; j++)
			{
				vold[i][j] = v[i][j];
			}
		}


		//////////////start solve continuity equation///////////
		//////////////start solve continuity equation///////////
		//////////////start solve continuity equation///////////
		//////////////start solve continuity equation///////////


		////////////////from left to right//////////////////////
		///////////////left side//////////////////////////////
		//////////////-as*Ts+ap/w*Tp-an*Tn=aw*Tw+ae*Te+b+ap*Tp*(1/w-1) 
		j = 1;
		for (i = 2; i < IPCV; i++)
		{
			dx1[i - 1] = 0;         //w
			dx2[i - 1] = rho * du[i][j] * dely;			//e
			ax[i - 1] = -rho * dv[i - 1][j] * delx;							//s
			cx[i - 1] = -rho * dv[i][j] * delx;				//n
			bx[i - 1] = -ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1];
			dx[i - 1] = dx2[i - 1] * pc[i][j + 1] + rho * (u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx;
		}

		i = 1;
		dx1[i - 1] = 0;         //w
		dx2[i - 1] = rho * du[i][j] * dely;			//e
		ax[i - 1] = 0;							//s
		cx[i - 1] = -rho * dv[i][j] * delx;				//n
		bx[i - 1] = -ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1];
		dx[i - 1] = dx2[i - 1] * pc[i][j + 1] + rho * (u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx;

		i = IPCV;
		dx1[i - 1] = 0;         //w
		dx2[i - 1] = rho * du[i][j] * dely;			//e
		ax[i - 1] = -rho * dv[i - 1][j] * delx;							//s
		cx[i - 1] = 0;				//n
		bx[i - 1] = -ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1];
		dx[i - 1] = dx2[i - 1] * pc[i][j + 1] + rho * (u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx;
		//					double *y = new double[N];
		TDMAI(y, ax, bx, cx, dx);
		for (k = 0; k < IPCV; k++)
		{
			pc[k + 1][j] = y[k];
		}

		/////////////////////////////middle ones///////////////////////////////////////

		for (j = 2; j < JPCV; j++)
		{
			for (i = 2; i < IPCV; i++)
			{
				dx1[i - 1] = rho * du[i][j - 1] * dely;         //w
				dx2[i - 1] = rho * du[i][j] * dely;			//e
				ax[i - 1] = -rho * dv[i - 1][j] * delx;							//s
				cx[i - 1] = -rho * dv[i][j] * delx;				//n
				bx[i - 1] = -ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1];
				dx[i - 1] = dx1[i - 1] * pc[i][j - 1] + dx2[i - 1] * pc[i][j + 1] + rho * (u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx;
			}

			i = 1;
			dx1[i - 1] = rho * du[i][j - 1] * dely;         //w
			dx2[i - 1] = rho * du[i][j] * dely;			//e
			ax[i - 1] = 0;							//s
			cx[i - 1] = -rho * dv[i][j] * delx;				//n
			bx[i - 1] = -ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1];
			dx[i - 1] = dx1[i - 1] * pc[i][j - 1] + dx2[i - 1] * pc[i][j + 1] + rho * (u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx;

			i = IPCV;
			dx1[i - 1] = rho * du[i][j - 1] * dely;         //w
			dx2[i - 1] = rho * du[i][j] * dely;			//e
			ax[i - 1] = -rho * dv[i - 1][j] * delx;							//s
			cx[i - 1] = 0;				//n
			bx[i - 1] = -ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1];
			dx[i - 1] = dx1[i - 1] * pc[i][j - 1] + dx2[i - 1] * pc[i][j + 1] + rho * (u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx;

			TDMAI(y, ax, bx, cx, dx);
			for (k = 0; k < IPCV; k++)
			{
				pc[k + 1][j] = y[k];
			}
		}
		///////////////////////right ones///////////////////////////
		j = JPCV;
		for (i = 2; i < IPCV; i++)
		{
			dx1[i - 1] = rho * du[i][j - 1] * dely;         //w
			dx2[i - 1] = 0;			//e
			ax[i - 1] = -rho * dv[i - 1][j] * delx;							//s
			cx[i - 1] = -rho * dv[i][j] * delx;				//n
			bx[i - 1] = -ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1];
			dx[i - 1] = dx1[i - 1] * pc[i][j - 1] + rho * (u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx;
		}

		i = 1;
		dx1[i - 1] = rho * du[i][j - 1] * dely;         //w
		dx2[i - 1] = 0;			//e
		ax[i - 1] = 0;							//s
		cx[i - 1] = -rho * dv[i][j] * delx;			//n
		bx[i - 1] = -ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1];
		dx[i - 1] = dx1[i - 1] * pc[i][j - 1] + rho * (u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx;

		i = IPCV;
		dx1[i - 1] = rho * du[i][j - 1] * dely;         //w
		dx2[i - 1] = 0;			//e
		ax[i - 1] = -rho * dv[i - 1][j] * delx;							//s
		cx[i - 1] = 0;				//n
		bx[i - 1] = -ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1];
		dx[i - 1] = dx1[i - 1] * pc[i][j - 1] + rho * (u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx;

		TDMAI(y, ax, bx, cx, dx);
		for (k = 0; k < IPCV; k++)
		{
			pc[k + 1][j] = y[k];
		}






		///////////now start correct p//////////////
		for (i = 1; i < IPCV + 1; i++)
		{
			for (j = 1; j < JPCV + 1; j++)
			{
				p[i][j] = p[i][j] + wp * pc[i][j];
			}
		}
		////////////correct u and v//////////////
		////du 76, dv 67/////////
		for (i = 1; i < IPCV + 1; i++)
		{
			for (j = 1; j < JPCV; j++)
			{
				u[i][j] = u[i][j] + du[i][j] * (pc[i][j] - pc[i][j + 1]);
			}
		}

		for (i = 1; i < IPCV; i++)
		{
			for (j = 1; j < JPCV + 1; j++)
			{
				v[i][j] = v[i][j] + dv[i][j] * (pc[i][j] - pc[i + 1][j]);
			}
		}

		/////////control p convergence//////
		for (i = 1; i < IPCV + 1; i++)
		{
			for (j = 1; j < JPCV + 1; j++)
			{
				Rp = Rp + fabs(rho*(u[i][j - 1] - u[i][j])*dely + rho * (v[i - 1][j] - v[i][j])*delx);
			}
		}
		Rp = Rp / rho / Vinf / 0.04;
		/////////control p convergence//////

		/////////control u convergence//////
		temp = 0;
		for (i = 1; i < IPCV + 1; i++)
		{
			for (j = 1; j < JPCV; j++)
			{
				Ru = Ru + (du[i][j] * (dely))*fabs(u[i][j] - uold[i][j]);   //where du[i][j]*(dely) equals to ap
				temp = temp + fabs(du[i][j] * (dely)*u[i][j]);
			}
		}
		Ru = Ru / temp;
		/////////control v convergence//////
		temp = 0;
		for (i = 1; i < IPCV; i++)
		{
			for (j = 1; j < JPCV + 1; j++)
			{
				Rv = Rv + (dv[i][j] * (delx))*fabs(v[i][j] - vold[i][j]);   //where dv[i][j]*(delx) equals to ap
				temp = temp + fabs(dv[i][j] * (delx)*v[i][j]);
			}
		}
		Rv = Rv / temp;

		//cout<<"overall-conv="<<Rp<<endl;
		cout << "v-conv=" << Rv << endl;

		in << itall << "," << log10(Rp) << "," << log10(Ru) << "," << log10(Rv) << ",";
		in << "\n";

	} while (Rp > 1e-6 || Ru > 1e-6 || Rv > 1e-6);
	////////////whole iteration end here//
	////////////whole iteration end here//
	////////////whole iteration end here//

	in.close();//close the file



	in; in.open("u5.csv", ios::ate);
	double ipt = 0;
	double Rt = 0;
	do
	{
		for (i = 0; i < IPCV + 2; i++)
		{
			for (j = 0; j < JPCV + 2; j++)
			{
				Told[i][j] = T[i][j];
			}
		}
		ipt = ipt + 1;
		Rt = 0;

		//left to right except outflow
			//////////////-as*Ts+ap/w*Tp-an*Tn=aw*Tw+ae*Te+b+ap*Tp*(1/w-1) 
		for (j = 1; j < JPCV; j++)
		{

			for (i = 1; i < IPCV + 1; i++)
			{
				ftw[i - 1][j - 1] = rho * cp*u[i][j - 1] * dely;
				Ptw[i - 1][j - 1] = ftw[i - 1][j - 1] / (K*dtw[i - 1][j - 1]);
				APtw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ptw[i - 1][j - 1]), 5));
				dx1[i - 1] = (K*dtw[i - 1][j - 1])*APtw[i - 1][j - 1] + max(ftw[i - 1][j - 1], 0);


				fte[i - 1][j - 1] = rho * cp*u[i][j] * dely;
				Pte[i - 1][j - 1] = fte[i - 1][j - 1] / (K*dte[i - 1][j - 1]);
				APte[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pte[i - 1][j - 1]), 5));
				dx2[i - 1] = (K*dte[i - 1][j - 1])*APte[i - 1][j - 1] + max(-fte[i - 1][j - 1], 0);

				fts[i - 1][j - 1] = rho * cp*v[i - 1][j] * delx;
				Pts[i - 1][j - 1] = fts[i - 1][j - 1] / (K*dts[i - 1][j - 1]);
				APts[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pts[i - 1][j - 1]), 5));
				ax[i - 1] = (K*dts[i - 1][j - 1])*APts[i - 1][j - 1] + max(fts[i - 1][j - 1], 0);
				ax[i - 1] = -ax[i - 1];

				ftn[i - 1][j - 1] = rho * cp*v[i][j] * delx;
				Ptn[i - 1][j - 1] = ftn[i - 1][j - 1] / (K*dtn[i - 1][j - 1]);
				APtn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ptn[i - 1][j - 1]), 5));
				cx[i - 1] = (K*dtn[i - 1][j - 1])*APtn[i - 1][j - 1] + max(-ftn[i - 1][j - 1], 0);
				cx[i - 1] = -cx[i - 1];

				bx[i - 1] = (-ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1]);
				dx[i - 1] = dx1[i - 1] * T[i][j - 1] + dx2[i - 1] * T[i][j + 1];
			}
			i = 1;
			dx[i - 1] = dx1[i - 1] * T[i][j - 1] + dx2[i - 1] * T[i][j + 1] + (-ax[i - 1])*T[i - 1][j];
			i = IPCV;
			dx[i - 1] = dx1[i - 1] * T[i][j - 1] + dx2[i - 1] * T[i][j + 1] + (-cx[i - 1])*T[i + 1][j];

			//////***********coefficient end******begin solve///////	
			double *y = new double[IPCV];
			TDMAI(y, ax, bx, cx, dx);
			for (k = 0; k < IPCV; k++)
			{
				T[k + 1][j] = y[k];
			}
			//delete[] y;


			for (i = 1; i < IPCV + 1; i++)
			{
				dT[i - 1][j - 1] = bx[i - 1];
			}


		}



		//////***********coefficient end******begin solve///////

		///outflow
		j = JPCV;
		for (i = 1; i < IPCV + 1; i++)
		{
			ftw[i - 1][j - 1] = rho * cp*u[i][j - 1] * dely;
			Ptw[i - 1][j - 1] = ftw[i - 1][j - 1] / (K*dtw[i - 1][j - 1]);
			APtw[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ptw[i - 1][j - 1]), 5));
			dx1[i - 1] = (K*dtw[i - 1][j - 1])*APtw[i - 1][j - 1] + max(ftw[i - 1][j - 1], 0);


			fte[i - 1][j - 1] = rho * cp*u[i][j] * dely;
			Pte[i - 1][j - 1] = fte[i - 1][j - 1] / (K*dte[i - 1][j - 1]);
			APte[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pte[i - 1][j - 1]), 5));
			dx2[i - 1] = (K*dte[i - 1][j - 1])*APte[i - 1][j - 1] + max(-fte[i - 1][j - 1], 0);
			dx2[i - 1] = 0;


			fts[i - 1][j - 1] = rho * cp*v[i - 1][j] * delx;
			Pts[i - 1][j - 1] = fts[i - 1][j - 1] / (K*dts[i - 1][j - 1]);
			APts[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Pts[i - 1][j - 1]), 5));
			ax[i - 1] = (K*dts[i - 1][j - 1])*APts[i - 1][j - 1] + max(fts[i - 1][j - 1], 0);
			ax[i - 1] = -ax[i - 1];

			ftn[i - 1][j - 1] = rho * cp*v[i][j] * delx;
			Ptn[i - 1][j - 1] = ftn[i - 1][j - 1] / (K*dtn[i - 1][j - 1]);
			APtn[i - 1][j - 1] = max(0, pow(1 - 0.1*fabs(Ptn[i - 1][j - 1]), 5));
			cx[i - 1] = (K*dtn[i - 1][j - 1])*APtn[i - 1][j - 1] + max(-ftn[i - 1][j - 1], 0);
			cx[i - 1] = -cx[i - 1];

			bx[i - 1] = (-ax[i - 1] - cx[i - 1] + dx1[i - 1] + dx2[i - 1]);
			dx[i - 1] = dx1[i - 1] * T[i][j - 1] + dx2[i - 1] * T[i][j + 1];
		}
		i = 1;
		dx[i - 1] = dx1[i - 1] * T[i][j - 1] + (-ax[i - 1])*T[i - 1][j];
		i = IPCV;
		dx[i - 1] = dx1[i - 1] * T[i][j - 1] + (-cx[i - 1])*T[i + 1][j];
		////////**********************///////////////
		for (i = 1; i < IPCV + 1; i++)
		{
			dT[i - 1][j - 1] = bx[i - 1];
		}
		////////***********************///////////////
		//////***********coefficient end******begin solve///////	
		double *y = new double[IPCV];
		TDMAI(y, ax, bx, cx, dx);
		for (k = 0; k < IPCV; k++)
		{
			T[k + 1][j] = y[k];
		}
		delete[] y;
		//////solve T end////////////////////////////
		j = JPCV + 1;
		for (i = 0; i < IPCV + 1; i++)
		{
			T[i][j] = T[i][j - 1];
		}


		for (i = 0; i < IPCV + 1; i++)
		{
			for (j = 0; j < JPCV + 1; j++)
			{
				Rt = Rt + fabs((T[i][j] - Told[i][j]) / T[i][j]);
			}
		}


		in << ipt << "," << Rt << ",";
		in << "\n";

	} while (Rt > 1e-6);
	cout << "Rt=" << Rt << endl;
	/////////////////////////////////////////////////////
	in.close();//close the file
	//////////////////////////////////////////////
	in; in.open("u6.csv", ios::ate);
	temp = 0;
	double uavg = 0;
	double h = 0;
	double Nu;
	j = 0;
	temp = 0;
	for (i = 1; i < IPCV + 1; i++)
	{
		temp = temp + u[i][j - 1] * T[i][j] * dely;
	}
	temp = temp / 0.02;
	uavg = 0;
	for (i = 1; i < IPCV + 1; i++)
	{
		uavg = uavg + u[i][j - 1];
	}
	uavg = uavg / IPCV;
	temp = temp / uavg;        //tmean
	h = 0;
	h = K * (T[IPCV][j] - 100) / (dely / 2) / (temp - 100);
	Nu = h * 2 * 0.02 / K;
	in << Nu << ",";


	for (j = 1; j < JPCV + 1; j++)
	{
		temp = 0;
		for (i = 1; i < IPCV + 1; i++)
		{
			temp = temp + (u[i][j] + u[i][j - 1]) / 2 * T[i][j] * dely;
		}
		temp = temp / 0.02;
		uavg = 0;
		for (i = 1; i < IPCV + 1; i++)
		{
			uavg = uavg + (u[i][j] + u[i][j - 1]) / 2;
		}
		uavg = uavg / IPCV;
		temp = temp / uavg;        //tmean
		h = 0;
		h = K * (T[IPCV][j] - 100) / (dely / 2) / (temp - 100);
		Nu = h * 2 * 0.02 / K;
		in << Nu << ",";
	}
	/////////////******************
	j = JPCV + 1;
	temp = 0;
	for (i = 1; i < IPCV + 1; i++)
	{
		temp = temp + u[i][j - 1] * T[i][j] * dely;
	}
	temp = temp / 0.02;
	uavg = 0;
	for (i = 1; i < IPCV + 1; i++)
	{
		uavg = uavg + u[i][j];
	}
	uavg = uavg / IPCV;
	temp = temp / uavg;        //tmean
	h = 0;
	h = K * (T[IPCV][j] - 100) / (dely / 2) / (temp - 100);
	Nu = h * 2 * 0.02 / K;
	in << Nu << ",";
	in.close();//close the file
///////////*/////////////////////////////////////////

///////////*/////////////////////////////////////////
	cout << "overall-conv=" << Rp << endl;
	cout << "u-conv=" << Ru << endl;
	cout << "v-conv=" << Rv << endl;
	in; in.open("u1.csv", ios::ate);
	i = IPCV / 2 + 1;
	for (j = 0; j < JPCV + 1; j++)
	{
		in << u[i][j] / Vinf << ",";         //first one
	}

	in.close();//close the file

	in; in.open("u2.csv", ios::ate);
	j = 24;
	in << "u" << ",";
	for (i = 0; i < IPCV + 2; i++)
	{
		in << u[i][j] << ",";         //second one
	}
	in << "\n";
	in << "v" << ",";
	for (i = 0; i < IPCV + 1; i++)
	{
		in << (v[i][j] + v[i][j - 1]) / 2 << ",";         //second one
	}
	in << "\n";
	in << "T" << ",";
	for (i = 0; i < IPCV + 2; i++)
	{
		in << (T[i][j] + T[i][j - 1]) / 2 << ",";         //second one
	}
	in.close();//close the file

	in; in.open("u3.csv", ios::ate);
	for (j = 1; j < JPCV + 1; j++)
	{
		temp = 0;
		for (i = 1; i < IPCV + 1; i++)
		{
			temp = temp + (u[i][j] + u[i][j - 1]) / 2 * T[i][j] * dely;
		}
		temp = temp / 0.02;
		uavg = 0;
		for (i = 1; i < IPCV + 1; i++)
		{
			uavg = uavg + (u[i][j] + u[i][j - 1]) / 2;
		}
		uavg = uavg / IPCV;
		temp = temp / uavg;        //tmean
		in << temp << ",";
	}
	in.close();//close the file
///////////////////////////////////
}



