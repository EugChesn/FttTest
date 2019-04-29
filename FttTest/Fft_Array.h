#pragma once
#include <vector>
using namespace std;
template<typename T>
class Fft_Array
{
public:
	T ** array;
	int n;
	int m;
public:
	Fft_Array(int n, int m)
	{
		this->n = n;
		this->m = m;
		array = new T*[n];
		array[0] = new T[m * n];
		for (int i = 1; i < n; i++)
			array[i] = array[0] + i * m;
	}
	Fft_Array(const vector<vector<T>> & vec)
	{
		this->n = vec.size();
		this->m = vec[0].size();
		array = new T*[n];
		array[0] = new T[m * n];
		for (int i = 1; i < n; i++)
			array[i] = array[0] + i * m;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
				array[i][j] = vec[i][j];
		}
	}
	Fft_Array(const Fft_Array & old_arr ,int new_n, int new_m)
	{
		this->n = new_n;
		this->m = new_m;
		array = new T*[n];
		array[0] = new T[m * n];
		for (int i = 1; i < n; i++)
			array[i] = array[0] + i * m;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				if (i < old_arr.n && j < old_arr.m)
					array[i][j] = old_arr.array[i][j];
				else
					array[i][j] = 0;
			}
		}
	}
	Fft_Array(const Fft_Array & img, const Fft_Array & filter) //multiply
	{
		this->n = img.n;
		this->m = img.m;
		array = new T*[n];
		array[0] = new T[m * n];
		for (int i = 1; i < n; i++)
			array[i] = array[0] + i * m;

		double normalize = n * m; //сразу нормализуем
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
				array[i][j] = img.array[i][j] * filter.array[i][j] / normalize;
		}
	}
	void rand_arr()
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
				array[i][j] = 1 + rand() % 20;
		}
	}
	void show()
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
				std::cout << array[i][j] << " ";

			std::cout << std::endl;
		}
		cout << endl;
	}
	void rot180() 
	{
		int center = n / 2;
		int center_j = m / 2;

		///////////////////////////////// Reverse rows
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < center_j; j++)
			{
				T b = array[i][j];
				array[i][j] = array[i][m - j - 1]; //n - i - 1
				array[i][m - j - 1] = b; //n - i - 1
			}
		}
		/////////////////////////////// Rotate 180
		for (int i = 0; i < center; i++)
		{
			for (int j = 0; j < m; j++)
			{
				T b = conj(array[i][j]);
				array[i][j] = conj(array[n - i - 1][j]);
				array[n - i - 1][j] = b;
			}
		}
	}
	~Fft_Array()
	{
		delete[] array[0];
		delete array;
	}
};

void fft(complex<double> *in, complex<double> * out, int H, int W)
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

void fft_double_to_complex(double * in,complex<double>*out,int H,int W)
{
	fftw_plan plan = fftw_plan_dft_r2c_2d(H, W,in, reinterpret_cast<fftw_complex*>(out), FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}
void ifft_complex_to_double(complex<double> * in, double *out, int H, int W)
{
	fftw_plan plan = fftw_plan_dft_c2r_2d(H, W, reinterpret_cast<fftw_complex*>(in), (out), FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}

