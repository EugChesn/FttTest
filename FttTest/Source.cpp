#include<iostream>
#include<ctime>
#include<complex>
#include<fftw3.h>

#define N 4096
#define M 3072 //размер картинки

#define N_f 315 //размер фильтра
#define M_f 315

using namespace std;
template<typename T>
void fft(T *in, T *out, int H, int W)
{
	fftw_plan plan = fftw_plan_dft_2d(H, W, reinterpret_cast<fftw_complex*> (in), reinterpret_cast<fftw_complex*> (in), FFTW_FORWARD, FFTW_ESTIMATE);
	//fftw_plan plan = fftw_plan_dft_c2r_2d(H, W, reinterpret_cast<double*> (in), reinterpret_cast<fftw_complex*> (out), FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}
void ifft(complex<double> *out, complex<double> *in, int H, int W)
{
	fftw_plan p = fftw_plan_dft_2d(H, W, reinterpret_cast<fftw_complex*> (in), reinterpret_cast<fftw_complex*> (in), FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}
template<typename T>
void show_arr(T ** arr, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			cout << arr[i][j] << " ";

		cout << endl;
	}
}
complex<double> ** create_array(int n,int m) 
{
	complex<double>**arr;
	arr = new complex<double>*[n];
	arr[0] = new complex<double>[m * n];
	for (int i = 1; i < n; i++)
		arr[i] = arr[0] + i * m;

	return arr;
}
void rand_arr(complex<double> ** arr , int H,int W)
{
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
			arr[i][j] = complex<double>(1 + rand() % 20, 1 + rand() % 20);
	}
}
void filling_zer(complex<double> ** new_arr, int new_size_n , int new_size_m , complex<double> ** old_arr ,int size_n ,int size_m) // заполнение нулями
{
	for (int i = 0; i < new_size_n; i++)
	{
		for (int j = 0; j < new_size_m; j++)
		{
			if (i < size_n && j < size_m)
				new_arr[i][j] = old_arr[i][j];
			else
				new_arr[i][j] = 0;
		}
	}
}
void normalize(complex<double> **in ,int H,int W)
{
	double norm = (H * W);
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
			in[i][j]/=norm;
	}
}
void multiply(complex<double> ** arr_ftt1, complex<double> ** arr_ftt2 ,complex<double> ** res , int H,int W)
{
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
			res[i][j] = arr_ftt1[i][j] * arr_ftt2[i][j];
	}
}
void rot180(complex<double> **arr, int n, int m) //for filter
{
	int center = n / 2; 
	int center_j = m / 2;

	///////////////////////////////// Reverse rows
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < center_j; j++) 
		{
			complex<double> b = arr[i][j];
			arr[i][j] = arr[i][m - j - 1]; //n - i - 1
			arr[i][m - j - 1] = b; //n - i - 1
		}
	}
	/////////////////////////////// Rotate 180
	for (int i = 0; i < center; i++)
	{
		for (int j = 0; j < m; j++)
		{
			complex<double> b = conj(arr[i][j]);
			arr[i][j] = conj(arr[n - i - 1][j]);
			arr[n - i - 1][j] = b; 
		}
	}
}
void delete_array(complex<double>** arr)
{
	delete[] arr[0];
	delete arr;
}
complex<double> ** xcorr2(complex<double> ** image, complex<double> ** filter ,int &n,int &m)
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

	complex<double>** mult = create_array(res_size_n, res_size_m); //Поэлементное перемножение результатов Фурье
	multiply(image_f, filter_f,mult, res_size_n, res_size_m); //out = in

	delete_array(image_f);
	delete_array(filter_f);

	ifft(&mult[0][0], &mult[0][0], res_size_n, res_size_m);//out = in //Обратное преобразование Фурье и нормализиция
	//normalize
	/*cout << "Result" << endl;
	show_arr(mult, res_size_n, res_size_m);*/
	return mult;
}

int main()
{
	srand(time(0));
	complex<double> **image = nullptr;
	complex<double> **filter = nullptr;
	image = create_array(N,M);
	filter = create_array(N_f,M_f);
	rand_arr(image,N,M);
	rand_arr(filter, N_f, M_f);

	
	//cout << "image" << endl;
	//show_arr(image,N,M);
	//cout << "filter" << endl;
	//show_arr(filter, N_f, M_f);

	int n_res = 0;
	int m_res = 0;
	complex<double>** res = xcorr2(image, filter,n_res,m_res);
	//show_arr(image, N, M);
	//cout << "end" << endl;

	//Test Furie
	/*fft(&image[0][0], &image[0][0], N, M);
	ifft(&image[0][0], &image[0][0], N, M);
	normalize(image, N, M);*/

	delete_array(res);

	delete_array(image);
	delete_array(filter);
	system("pause");
	return 0;
}