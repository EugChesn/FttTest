#include<iostream>
#include<ctime>
#include<complex>
#include<fftw3.h>
#include"Fft_Array.h"

#define N 3
#define M 3 //размер картинки

#define N_f 3 //размер фильтра
#define M_f 3

using namespace std;
/*complex<double> ** xcorr2(complex<double> ** image, complex<double> ** filter ,int &n,int &m)
{
	int res_size_n = N + N_f - 1; n = res_size_n;
	int res_size_m = M + M_f - 1; m = res_size_m;
	complex<double>** image_f = create_array(res_size_n, res_size_m);
	complex<double>** filter_f = create_array(res_size_n, res_size_m);
	rot180(filter,N_f,M_f);

	filling_zer(image_f, res_size_n, res_size_m, image, N, M); // Дополняем нулями
	filling_zer(filter_f, res_size_n, res_size_m, filter, N_f, M_f);

	fft(&image_f[0][0], &image_f[0][0], res_size_n, res_size_m); //out = in //Фурье прямое преобразование
	fft(&filter_f[0][0], &filter_f[0][0], res_size_n, res_size_m); //out = in

	complex<double>** mult = create_array(res_size_n, res_size_m); //Поэлементное перемножение результатов Фурье и нормализиция
	multiply(image_f, filter_f,mult, res_size_n, res_size_m); //out = in

	delete_array(image_f);
	delete_array(filter_f);

	ifft(&mult[0][0], &mult[0][0], res_size_n, res_size_m);//out = in //Обратное преобразование Фурье
	return mult;
}*/
template<typename T>
Fft_Array<T> * xcorr2(Fft_Array<T> & image, Fft_Array<T> & filter)
{
	filter.rot180();
	int res_size_n = image.n + filter.n - 1; 
	int res_size_m = image.m + filter.n - 1; 

	Fft_Array<T> image_temp(image, res_size_n, res_size_m);
	Fft_Array<T> filter_temp(filter, res_size_n, res_size_m);

	fft(&image_temp.array[0][0], &image_temp.array[0][0], res_size_n, res_size_m); //out = in //Фурье прямое преобразование
	fft(&filter_temp.array[0][0], &filter_temp.array[0][0], res_size_n, res_size_m); //out = in

	Fft_Array<T> * mult = new Fft_Array<T>(image_temp, filter_temp);
	ifft(&mult->array[0][0], &mult->array[0][0], res_size_n, res_size_m);//out = in //Обратное преобразование Фурье

	return mult;
}
int main()
{
	srand(time(0));

	//Xcorr2 for complex to complex 
	Fft_Array<complex<double>> image(N_f, M_f);
	Fft_Array<complex<double>> filter(N, M);
	image.rand_arr();
	filter.rand_arr();
	Fft_Array<complex<double>> *result = xcorr2(image, filter);


	/*Fft_Array<double> image(N, M);
	Fft_Array<double> filter(N_f, M_f);
	image.rand_arr();
	filter.rand_arr();
	image.rand_arr();
	filter.rand_arr();
	image.show();
	filter.show();
	int size_complex = floor(M / 2) + 1;
	Fft_Array<complex<double>> tmp(N, size_complex);
	//fftw_plan plan = fftw_plan_dft_r2c_2d(N ,M, &image.array[0][0], reinterpret_cast<fftw_complex*>(&tmp.array[0][0]), FFTW_ESTIMATE);
	Fft_Array<double> out(N, M);
	Fft_Array<complex<double>> in(N, M);
	in.rand_arr();
	in.show();
	fftw_plan p = fftw_plan_dft_c2r_2d(N, M, reinterpret_cast<fftw_complex*> (&in.array[0][0]), (&out.array[0][0]), FFTW_ESTIMATE);
	fftw_execute(p);
	out.show();*/


	cout << "runtime = " << clock() / 1000.0 << endl; // время работы программы 
	system("pause");
	return 0;
}