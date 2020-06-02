
#include <stdio.h>
#include <math.h>
#include <iostream>

double Psubi(int i, int k, double* e, double* mu, int* r);

double Psmoi(int i, int num_vymog, int cn, double* e, double* mu, int* r, int N);

double Li(int i, double cn, double* e, double* mu, int* r, int N);

double Ri(int i, double cn, double* e, double* mu, int* r, int N);

double fact(double x) {
	if (x == 0 || x == 1) {
		return 1;
	}
	else {
		return x*fact(x - 1);
	}
}

int main()
{
	const int N = 20;
	double e[] = { 1, 1.25 };
	double inverse_mu[] = { 0.8, 0.3};
	int r[] = { 1, 3 };
	double CN = 0;
	for (int num_vymog = 0; num_vymog <= N; num_vymog++) {
		CN += Psubi(0, num_vymog, e, inverse_mu, r) * Psubi(1, N - num_vymog, e, inverse_mu, r);
	}
	CN = pow(CN, -1);
	printf("Normalizing coef: %f;\n", CN);
	double test = 0;
	for (int num_vymog = 0; num_vymog <= N; num_vymog++) {
		test += Psmoi(0, num_vymog, CN, e, inverse_mu, r, N);
	}
	printf("Sum of probs (test model credibility): %f; \n", test);
	double L1 = Li(0, CN, e, inverse_mu, r, N);
	double L2 = Li(1, CN, e, inverse_mu, r, N);
	double R1 = Ri(0, CN, e, inverse_mu, r, N);
	double R2 = Ri(1, CN, e, inverse_mu, r, N);
	printf("L1: %f; \n", L1);
	printf("L2: %f; \n", L2);
	printf("R1: %f; \n", R1);
	printf("R2: %f; \n", R2);
	printf("M1: %f; \n", L1 + R1);
	printf("M2: %f; \n", L2 + R2);
	printf("Total expectation of requests in network: %f; \n", L1 + R1 + L2 + R2);
	printf("Intensity 1 (lambda_1): %f; \n", R1 / inverse_mu[0]);
	printf("Intensity 2 (lambda_2): %f; \n", R2 / inverse_mu[1]);
	printf("Mean time in system 1 (T_1): %f; \n", (L1 + R1) / (R1 / inverse_mu[0]));
	printf("Mean time in system 2 (T_2): %f; \n", (L2 + R2) / (R2 / inverse_mu[1]));
	printf("Mean time in queue in system 1 (T_1): %f; \n", L1 / (R1 / inverse_mu[0]));
	printf("Mean time in queue in system 2 (T_2): %f; \n", L2 / (R2 / inverse_mu[1]));
}


double Psubi(int smo_id, int num_vymog, double* e, double* mu, int* num_channels) {
	double second = 0;
	if (num_vymog > num_channels[smo_id]) {
		second = 1 / (fact(num_channels[smo_id]) * pow(num_channels[smo_id], num_vymog - num_channels[smo_id]));
	}
	else {
		second = 1 / fact(num_vymog);
	}
	return pow((e[smo_id] * mu[smo_id]) , num_vymog) * second;
}

double Psmoi(int i, int num_vymog, int cn, double* e, double* mu, int* r, int N) {
	int ii = 0;
	if (i == 0)
	{
		ii = 1;
	}
	return cn * Psubi(i, num_vymog, e, mu, r) * Psubi(ii, N-num_vymog, e, mu, r);
}

double Li(int i, double cn, double* e, double* mu, int* r, int N) {
	double mexp = 0;
	for (int j = r[i] + 1; j <= N; j++) {
		mexp += (j - r[i]) * Psmoi(i, j, cn, e, mu, r, N);
	}
	return mexp;
}

double Ri(int i, double cn, double* e, double* mu, int* r, int N) {
	double mexp = 0;
	for (int j = 0; j <= r[i]-1; j++) {
		mexp += (r[i] - j) * Psmoi(i, j, cn, e, mu, r, N);
	}
	return r[i] - mexp;
}