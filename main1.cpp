#include <stdio.h>
#pragma warning(disable:4996)
#include <string.h>
#include "matrix.h"
#include <stdlib.h>

matrix		invA;
matrix		invB;
matrix		D;
matrix		E;
double		sigma_1_2_r;
double		sigma_1_2_u;
double	D1_r_u;
double	D1_r_v;
double	D1_v_r;
double	D1_r;
double	D1_u;
double	D2_r_u;
double	D2_r_v;
double	D2_v_r;
double	D2_r;
double	D2_u;
double	h2;
double	sigma_a1_r;
double	sigma_a2_r;
double	sigma_a1_u;
double	sigma_a2_u;

int get_invA_invB()
{

	int		x, y;
	char*   	m_map;	// material map
	int*    	t_map;	// type map
	double* 	pA;
	double* 	pB;
	char*  	 	pS;
	int		N;
	int		N2;
	int		k;


	printf("N:");
	scanf("%d", &N);
	char reflector[10];
	printf("material (H2O,D2O,Be): ");
	scanf("%s", reflector);
	N2 = N * N;

	D1_r_u = 1.0;
	D1_r_v = 1.0;
	D1_v_r = 1.0;
	D1_r = 1.0;
	D1_u = 0.7896;
	D2_r_u = 1.0;
	D2_r_v = 1.0;
	D2_v_r = 1.0;
	D2_r = 1.0;
	D2_u = 0.7848;
	h2 = (0.2/(double)N)/(0.2/(double)N);


	sigma_a1_u = 4.398*0.01;
	sigma_1_2_u = 1.098*0.001;
	sigma_a2_u = 6*0.1;

	if (!strcmp(reflector, "H2O"))
	{
		D1_r_u =0.4324;
		D1_r_v = 0.0375;
		D1_v_r = 0.0375;
		D1_r = 0.0751;

		D2_r_u = 0.4044;
		D2_r_v = 0.0119;
		D2_v_r = 0.0119;
		D2_r = 0.0239;

		sigma_a1_r  = 6.523*0.001;
		sigma_1_2_r = 9.449*0.1;
		sigma_a2_r  = 2.521*0.1;
	}
	else if (!strcmp(reflector, "D2O"))
	{
		D1_r_u = 0.4764;
		D1_r_v = 0.0816;
		D1_v_r = 0.0816;
		D1_r = 0.1632;

		D2_r_u = 0.4416;
		D2_r_v = 0.0492;
		D2_v_r = 0.0492;
		D2_r = 0.0983;

		sigma_a1_r = 1.213*0.00001;
		sigma_1_2_r = 1.094*0.1;
		sigma_a2_r = 3.882*0.0001;
	}
	else if (!strcmp(reflector, "Be"))
	{
		D1_r_u = 0.4300;
		D1_r_v = 0.0327;
		D1_v_r = 0.0327;
		D1_r = 0.0705;

		D2_r_u =0.4206;
		D2_r_v =0.0282;
		D2_v_r = 0.0282;
		D2_r = 0.0563;

		sigma_a1_r = 6.623*0.001;
		sigma_1_2_r = 26.358*0.01;
		sigma_a2_r = 6.629*0.001;
	}
	else
	{
		D1_r_u = 0.4384;
		D1_r_v = 0.0436;
		D1_v_r = 0.0436;
		D1_r = 0.0871;

		D2_r_u = 0.4287;
		D2_r_v = 0.0145;
		D2_v_r = 0.0145;
		D2_r = 0.0726;

		sigma_a1_r = 2.330*0.0001;
		sigma_1_2_r = 3.443*0.01;
		sigma_a2_r = 2.594*0.001;
	}


	m_map = new char[N * N];
	t_map = new int[N * N];
	pA = new double[N2 * N2];
	pB = new double[N2 * N2];
	pS = new char[N2 * N2];

	for (y = 0; y < N; y++) {
		for (x = 0; x < N; x++) {
			m_map[y * N + x] = 1;
		}
	}

	for (y = 0; y < N; y++) {
		m_map[y * N + N - 1] = 2;
	}

	for (x = 0; x < N; x++) {
		m_map[y * (N-1) + x] = 2;
	}



	for (y = 0; y < N; y++)
		for (x = 0; x < N; x++) {
			if (y == 0) {
				if (x == 0) t_map[y * N + x] = 11;
				else if(x == N - 1) t_map[y * N + x] = 14;
				else if(x == N - 2) t_map[y * N + x] = 13;
				else			    t_map[y * N + x] = 12;
			}
			else if (y == N-1) {
				if (x == 0) t_map[y * N + x] = 1;
				else if (x == N - 1) t_map[y * N + x] = 3;
				else			     t_map[y * N + x] = 2;
			}
			else if (y == N-2) {
				if (x == 0) t_map[y * N + x] = 4;
				else if (x == N - 1) t_map[y * N + x] = 10;
				else if (x == N - 2) t_map[y * N + x] = 6;
				else			  t_map[y * N + x] = 5;
			}
			else {
				if (x == 0) t_map[y * N + x] = 7;
				else if (x == N - 1) t_map[y * N + x] = 10;
				else if (x == N - 2) t_map[y * N + x] = 9;
				else			  t_map[y * N + x] = 8;
			}
		}

	printf("------------material map---------------\n");
	printf("  1 : uranium\n");
	printf("  2 : reflector\n");

	for (y = 0; y < N; y++) {
		for (x = 0; x < N; x++) {
			printf("%02d ", m_map[y * N + x]);
		}
		printf("\n");
	}


	printf("------------coeffients type map---------------\n");
	printf("  case 1 ~ 14\n");

	for (y = 0; y < N; y++) {
		for (x = 0; x < N; x++) {
			printf("%02d ", t_map[y * N + x]);
		}
		printf("\n");
	}


	int			t_index;
	double		a_coeff_1[14];
	double		b_coeff_1[14];
	double		c_coeff_1[14];
	double		d_coeff_1[14];
	double		e_coeff_1[14];

	double		a_coeff_2[14];
	double		b_coeff_2[14];
	double		c_coeff_2[14];
	double		d_coeff_2[14];
	double		e_coeff_2[14];

	for (k = 0; k < 14; k++) {
		a_coeff_1[k] = 0;
		b_coeff_1[k] = 0;
		c_coeff_1[k] = 0;
		d_coeff_1[k] = 0;
		e_coeff_1[k] = 0;

		a_coeff_2[k] = 0;
		b_coeff_2[k] = 0;
		c_coeff_2[k] = 0;
		d_coeff_2[k] = 0;
		e_coeff_2[k] = 0;
	}


	D1_r_u *= 2.2;
	D1_r_v *= 2.2;
	D1_v_r *= 2.2;
	D1_r   *= 2.2;
	D1_u   *= 2.2;
	D2_r_u *= 2.2;
	D2_r_v *= 2.2;
	D2_v_r *= 2.2;
	D2_r   *= 2.2;
	D2_u   *= 2.2;

	sigma_a1_r *= 2.4;
	sigma_a2_r *= 2.4;
	sigma_1_2_r *= 2.4;
	sigma_a1_u *= 2.4;
	sigma_a2_u *= 2.4;
	sigma_1_2_u *= 2.4;


	//  case 1
	d_coeff_1[0] = -(D1_r_u / h2);
	b_coeff_1[0] = -(D1_r / h2) + (4 * D1_v_r / h2) + (D1_r_u / h2) + sigma_a1_r + sigma_1_2_r;
	c_coeff_1[0] = -(D1_r / h2);


	d_coeff_2[0] = -(D2_r_u / h2);
	b_coeff_2[0] = -(D2_r / h2) + (4 * D2_v_r / h2) + (D2_r_u / h2) + sigma_a2_r;
	c_coeff_2[0] = -(D2_r / h2);

	//  case 2
	d_coeff_1[1] = -(D1_r_u / h2);
	a_coeff_1[1] = -(D1_r / h2);
	b_coeff_1[1] = (4 * D1_r_v / h2) + (D1_r_u / h2) + (2 * D1_r / h2) + sigma_a1_r + sigma_1_2_r;
	c_coeff_1[1] = -(D1_r / h2);


	d_coeff_2[1] = -(D2_r_u / h2);
	a_coeff_2[1] = -(D2_r / h2);
	b_coeff_2[1] = (4 * D2_r_v / h2) + (D2_r_u / h2) + (2 * D2_r / h2) + sigma_a2_r;
	c_coeff_2[1] = -(D2_r / h2);


	//	case 3
	d_coeff_1[2] = -(D1_r / h2);
	a_coeff_1[2] = -(D1_r / h2);
	b_coeff_1[2] = (8 * D1_r_v / h2) + (2 * D1_r / h2) + sigma_a1_r + sigma_1_2_r;


	d_coeff_2[2] = -(D2_r / h2);
	a_coeff_2[2] = -(D2_r / h2);
	b_coeff_2[2] = (8 * D2_r_v / h2) + (2 * D2_r / h2) + sigma_a2_r;


	//	case 4
	d_coeff_1[3] = -(D1_u / h2);
	b_coeff_1[3] = (2 * D1_u / h2) + (D1_r_u / h2) + sigma_a1_u + sigma_1_2_u;
	c_coeff_1[3] = -(D1_u / h2);
	e_coeff_1[3] = -(D1_r_u / h2);


	d_coeff_2[3] = -(D2_u / h2);
	b_coeff_2[3] = (2 * D2_u / h2) + (D2_r_u / h2) + sigma_a2_u;
	c_coeff_2[3] = -(D2_u / h2);
	e_coeff_2[3] = -(D2_r_u / h2);

	//	case 5 
	d_coeff_1[4] = -(D1_u / h2);
	a_coeff_1[4] = -(D1_u / h2);
	b_coeff_1[4] = (3 * D1_u / h2) + (D1_r_u / h2) + sigma_a1_u + sigma_1_2_u;
	c_coeff_1[4] = -(D1_u / h2);
	e_coeff_1[4] = -(D1_r_u / h2);

	d_coeff_2[4] = -(D2_u / h2);
	a_coeff_2[4] = -(D2_u / h2);
	b_coeff_2[4] = (3 * D2_u / h2) + (D2_r_u / h2) + sigma_a2_u;
	c_coeff_2[4] = -(D2_u / h2);
	e_coeff_2[4] = -(D2_r_u / h2);

	//	case 6
	d_coeff_1[5] = -(D1_u / h2);
	a_coeff_1[5] = -(D1_u / h2);
	b_coeff_1[5] = (2 * D1_u / h2) + (2 * D1_r_u / h2) + sigma_a1_u + sigma_1_2_u;
	c_coeff_1[5] = -(D1_r_u / h2);
	e_coeff_1[5] = -(D1_r_u / h2);

	d_coeff_2[5] = -(D2_u / h2);
	a_coeff_2[5] = -(D2_u / h2);
	b_coeff_2[5] = (2 * D2_u / h2) + (2 * D2_r_u / h2) + sigma_a2_u;
	c_coeff_2[5] = -(D2_r_u / h2);
	e_coeff_2[5] = -(D2_r_u / h2);

	//	case 7 
	d_coeff_1[6] = -(D1_u / h2);
	b_coeff_1[6] = (3 * D1_u / h2) + sigma_a1_u + sigma_1_2_u;
	c_coeff_1[6] = -(D1_u / h2);
	e_coeff_1[6] = -(D1_u / h2);

	d_coeff_2[6] = -(D2_u / h2);
	b_coeff_2[6] = (3 * D2_u / h2) + sigma_a2_u;
	c_coeff_2[6] = -(D2_u / h2);
	e_coeff_2[6] = -(D2_u / h2);


	//	case 8
	d_coeff_1[7] = -(D1_u / h2);
	a_coeff_1[7] = -(D1_u / h2);
	b_coeff_1[7] = (4 * D1_u / h2) + sigma_a1_u + sigma_1_2_u;
	c_coeff_1[7] = -(D1_u / h2);
	e_coeff_1[7] = -(D1_u / h2);

	d_coeff_2[7] = -(D2_u / h2);
	a_coeff_2[7] = -(D2_u / h2);
	b_coeff_2[7] = (4 * D2_u / h2) + sigma_a2_u;
	c_coeff_2[7] = -(D2_u / h2);
	e_coeff_2[7] = -(D2_u / h2);


	//	case 9
	d_coeff_1[8] = -(D1_u / h2);
	a_coeff_1[8] = -(D1_u / h2);
	b_coeff_1[8] = (3 * D1_u / h2) + (D1_r_u / h2) + sigma_a1_u + sigma_1_2_u;
	c_coeff_1[8] = -(D1_r_u / h2);
	e_coeff_1[8] = -(D1_u / h2);


	d_coeff_2[8] = -(D2_u / h2);
	a_coeff_2[8] = -(D2_u / h2);
	b_coeff_2[8] = (3 * D2_u / h2) + (D2_r_u / h2) + sigma_a2_u;
	c_coeff_2[8] = -(D2_r_u / h2);
	e_coeff_2[8] = -(D2_u / h2);

	//	case 10 
	d_coeff_1[9] = -(D1_r / h2);
	a_coeff_1[9] = -(D1_r_u / h2);
	b_coeff_1[9] = (4 * D1_r_v / h2) + (D1_r_u / h2) + (2 * D1_r / h2) + sigma_a1_r + sigma_1_2_r;
	e_coeff_1[9] = -(D1_r / h2);

	d_coeff_2[9] = -(D2_r / h2);
	a_coeff_2[9] = -(D2_r_u / h2);
	b_coeff_2[9] = (4 * D2_r_v / h2) + (D2_r_u / h2) + (2 * D2_r / h2) + sigma_a2_r;
	e_coeff_2[9] = -(D2_r / h2);



	//	case 11
	b_coeff_1[10] = (2 * D1_u / h2) + sigma_a1_u + sigma_1_2_u;
	c_coeff_1[10] = -(D1_u / h2);
	e_coeff_1[10] = -(D1_u / h2);

	b_coeff_2[10] = (2 * D2_u / h2) + sigma_a2_u;
	c_coeff_2[10] = -(D2_u / h2);
	e_coeff_2[10] = -(D2_u / h2);

	//	case 12
	a_coeff_1[11] = -(D1_u / h2);
	b_coeff_1[11] = (3 * D1_u / h2) + sigma_a1_u + sigma_1_2_u;
	c_coeff_1[11] = -(D1_u / h2);
	e_coeff_1[11] = -(D1_u / h2);

	a_coeff_2[11] = -(D2_u / h2);
	b_coeff_2[11] = (3 * D2_u / h2) + sigma_a2_u;
	c_coeff_2[11] = -(D2_u / h2);
	e_coeff_2[11] = -(D2_u / h2);


	//	case 13
	a_coeff_1[12] = -(D1_u / h2);
	b_coeff_1[12] = (2 * D1_u / h2) + (D1_r_u / h2) + sigma_a1_u + sigma_1_2_u;
	c_coeff_1[12] = -(D1_r_u / h2);
	e_coeff_1[12] = -(D1_u / h2);

	a_coeff_2[12] = -(D2_u / h2);
	b_coeff_2[12] = (2 * D2_u / h2) + (D2_r_u / h2) + sigma_a2_u;
	c_coeff_2[12] = -(D2_r_u / h2);
	e_coeff_2[12] = -(D2_u / h2);

	//	case 14
	a_coeff_1[13] = -(D1_r_u / h2);
	b_coeff_1[13] = (8 * D1_r_v / h2) + (D1_r_u / h2) + (D1_r / h2) + sigma_a1_r + sigma_1_2_r;
	e_coeff_1[13] = -(D1_r / h2);

	a_coeff_2[13] = -(D2_r_u / h2);
	b_coeff_2[13] = (8 * D2_r_v / h2) + (D2_r_u / h2) + (D2_r / h2) + sigma_a2_r;
	e_coeff_2[13] = -(D2_r / h2);




	for (y = 0; y < N2; y++) {
		for (x = 0; x < N2; x++) {
			pA[y * N2 + x] = 0.0;
			pB[y * N2 + x] = 0.0;

			pS[y * N2 + x] = ' ';
		}
		t_index = t_map[y] - 1;



		pA[y * N2 + y] = b_coeff_1[t_index];
		pB[y * N2 + y] = b_coeff_2[t_index];
		pS[y * N2 + y] = 'b';
		if ((y - 1) >= 0) {
			pA[y * N2 + y - 1] = a_coeff_1[t_index];
			pB[y * N2 + y - 1] = a_coeff_2[t_index];
			pS[y * N2 + y - 1] = 'a';
		}
		if ((y + 1) < N * N) {
			pA[y * N2 + y + 1] = c_coeff_1[t_index];
			pB[y * N2 + y + 1] = c_coeff_2[t_index];
			pS[y * N2 + y + 1] = 'c';
		}
		if ((y + N) < N * N) {
			pA[y * N2 + y + N] = e_coeff_1[t_index];
			pB[y * N2 + y + N] = e_coeff_2[t_index];
			pS[y * N2 + y + N] = 'e';
		}
		if ((y - N) >= 0) {
			pA[y * N2 + y - N] = d_coeff_1[t_index];
			pB[y * N2 + y - N] = d_coeff_2[t_index];
			pS[y * N2 + y - N] = 'd';
		}
	}

	printf("----------coefficients map---------------\n");
	for (y = 0; y < N2; y++) {
		for (x = 0; x < N2; x++) {
			printf("%c ", pS[y * N2 + x]);
		}
		printf("  %d \n", t_map[y]);
	}

	matrix		A(N2, N2, pA);
	matrix		B(N2, N2, pB);

	printf("------------------A---------------------\n");
	A.print();

	invA.set_size(N2, N2);
	invB.set_size(N2, N2);

	D.set_size(N2, 1);
	E.set_size(N2, 1);

	for(y=0; y <N; y++)
		for(x=0; x<N; x++) {	
			if(m_map[y*N + x] == 1) {
				D.set_value(y*N + x, 0, sigma_1_2_u);
			} else {
				D.set_value(y*N + x, 0, sigma_1_2_r);
			}
		}



	matrix::inverse(A, invA);
	matrix::inverse(B, invB);

	return N;
}


int main()
{
	int			N;
	int			N2;	
	int			i;
	int			n;
	double		k;
	double		t;
	int			x, y;
	double		v_sigma_f1;
	double		v_sigma_f2;
	double		sum_p;
	double		sum_c;
	

	N = get_invA_invB();
	N2 = N*N;	

	v_sigma_f1 = 3.758*0.01;
	v_sigma_f2 = 1.158;
	k = 1.0;

	matrix		phi_1_i(N2, 1); 
	matrix		n_phi_1_i(N2, 1);
	matrix		phi_2_i(N2, 1);
	matrix		C(N2, 1);






	for(y=0;y<N; y++) {
		for(x=0; x<N; x++) {
			if((x < N-1) && (y < N-1)) {
				phi_1_i.set_value(y*N + x, 0, 1.0);
				phi_2_i.set_value(y*N + x, 0, 1.0);
			} else {
				phi_1_i.set_value(y*N + x, 0, 0.0);
				phi_2_i.set_value(y*N + x, 0, 0.0);
			}	 
/*
			if(((y == 0)  && (x < N-1)) ||
		       ((x == 0) && (y < N-1))) {
				phi_1_i.set_value(y*N + x, 0, 1.0/(2*N-1));
				phi_2_i.set_value(y*N + x, 0, 1.0/(2*N-1));
			} else {
				phi_1_i.set_value(y*N + x, 0, 0);
				phi_2_i.set_value(y*N + x, 0, 0);
			}
*/
		}
	}


	sum_p = 0.0;
	for(i=0; i<N2; i++) {
		sum_p = sum_p + (v_sigma_f1*phi_1_i.get_value(i, 0)) + (v_sigma_f2*phi_2_i.get_value(i, 0));
	}
	sum_c = sum_p;
	
	char		fname[4096];

	for(n=0; n<32; n++) {
		printf("sum_p = %lf  sum_c = %lf k = %lf\n", sum_p, sum_c, k);
//		printf("------------ph_1-------------\n");
//		phi_1_i.print2(N, N);
		memset(fname, 0, 4096);
		sprintf(fname, "phi_1_%02d.dat", n);
		phi_1_i.save2(N, N, fname);


//		printf("------------ph_2-------------\n");
//		phi_2_i.print2(N, N);
		sprintf(fname, "phi_2_%02d.dat", n);
		phi_2_i.save2(N, N, fname);


		sum_p = 0.0;
		for(i=0; i<N2; i++) {
			sum_p = sum_p + (v_sigma_f1*phi_1_i.get_value(i, 0)) + (v_sigma_f2*phi_2_i.get_value(i, 0));
		}

		C = ((v_sigma_f1/k)*phi_1_i) + ((v_sigma_f2/k)*phi_2_i);



		n_phi_1_i = invA*C;






		for(i=0; i<N2; i++) {
			E.set_value(i, 0, phi_1_i.get_value(i, 0)*D.get_value(i, 0));
		}
		

		phi_1_i = n_phi_1_i;
		phi_2_i = invB*E;


		sum_c = 0.0;
		for(i=0; i<N2; i++) {
			sum_c = sum_c + (v_sigma_f1*phi_1_i.get_value(i, 0)) + (v_sigma_f2*phi_2_i.get_value(i, 0));
		}

		
		k = k*(sum_c/sum_p);	// k(i) = k(i-1)*sum_c/sum_p;
	}
	
}
