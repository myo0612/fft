#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <math.h>
#define DATASIZE 64
#define PI 3.1415926536 // 円周率
#define DELTA 0.2 // サンプリング周期

// フーリエ変換
//******************************************************************************************************
// data[DATASIZE]を高速フーリエ変換した結果を引数data[DATASIZE]に格納し直す
// DATASIZEとDELTA(dt)は固定なので外から入力できるよう対応する必要あり
// 下記URLのExcelマクロを参考にC言語へ書き直した。Excelにて使用実績あり、出力値が一致したことを確認済。
//		http://tsuyu.cocolog-nifty.com/blog/2007/03/publi.html
// 処理内容：
// 				入力データを原点へ合わせる
// 				FFT
// 				FFT結果(複素数)の絶対値を取り、入力データ値の大きさと合わせたものを出力
//
//******************************************************************************************************

void dft(double data[DATASIZE] /*[入力] 実数数値列としての信号: 配列)*/ ){
	// 配列の要素数を格納
	// const int N = sizeof data /sizeof data[0];
	int N = DATASIZE;
	int N_half = N / 2;
	int k, p = 0;
	int m = log(N) / log(2); 
	double sum = 0;
	double avg = 0;
	double xr[N], xi[N], xd, s[N_half], c[N_half];

	double a = 0;
	double b = PI * 2 / N;

	// initializing data
	for (int i = 0; i < N; i++) {
		sum += data[i];
	}
	avg = sum / N;
	for (int i = 0; i < N; i++) {
		xr[i] = data[i] - avg;
		xi[i] = 0;
	}

	// creating array for FFT 
	for(int n = 0; n < N / 2; n++) {
		s[n] = sin(a);
		c[n] = cos(a);
		a    = a + b;
	}

	//=================================================
	// START FFT
	//=================================================
	int l = N;
	int h = 1;
	int j = 0;

	for (int g = 1; g <= m; g++) {
		l = l / 2;
		k = 0;

		for (int q = 1; q <= h; q++ ) {
			p = 0;

			for (int i = k; i < l + k; i++) {
				j = i + l;
				a = xr[i] - xr[j];
				b = xi[i] - xi[j];
				xr[i] = xr[i] + xr[j];
				xi[i] = xi[i] + xi[j];

				if (p == 0) {
					xr[j] = a;
					xi[j] = b;
				} else {
					xr[j] = a * c[p] + b * s[p];
					xi[j] = b * c[p] - a * s[p];
				}

				p = p + h;
			} // end of loop "i"

			k += 2 * l;
		} // end of loop "q"

		h += h;
	} // end of loop "g"

	j = N / 2;

	for (int i = 1; i < N - 1; i++) {
		k = N;
		if (j < i) {
			xd = xr[i];
			xr[i] = xr[j];
			xr[j] = xd;
			xd = xi[i];
			xi[i] = xi[j];
			xi[j] = xd;
		}
		k = k / 2;
		while (j >= k) {
			j = j - k;
			k = k / 2;
		}
		j = j + k;
	} // end of loop "i"

	// output data, IMABS/(data num)/2
	for (int i = 0; i < N / 2; i++) {
		data[i] = sqrt( pow(xr[i], 2) + pow(xi[i], 2) );
		data[i] = data[i] / (N / 2);
	}

} // end of function

// ディレクトリ内からファイルを開く
int main(int argc, char *argv[])
{

	DIR *dir;
	FILE *fp;
	struct dirent *dp;
	char path[128] = "/HOME/fft/data_dir/";
	char full_path[1024] = "";

	if ((dir = opendir(path)) == NULL) {
		fprintf(stderr, "Can't open directory %s\n", path);
		return -1;
	}

	int file_count = 0;

	for (dp = readdir(dir); dp != NULL; dp = readdir(dir)) {

		if ((strcmp(dp->d_name, ".") == 0) || (strcmp(dp->d_name, "..") == 0)) continue; 
		sprintf(full_path, "%s%s", path, dp->d_name);
		printf("Open a file: %s\n", full_path);
		if ((fp = fopen(full_path, "r")) == NULL) {
			fprintf(stderr, "%sのオープンに失敗しました.\n", full_path);
		}	
		int i=0;
		double data[DATASIZE];
		/* ファイルの終端まで文字を読み取り表示する */
		// while ( fgets(readline, DATASIZE, fp) != NULL ) {
		while( fscanf(fp,"%lf",&data[i]) != EOF ){
			// printf("data[%d] = %f\n", i, data[i]);
			i++;
		}
		FILE *wfp;
		char wfilename[1024] = "";
		char *base_wfilename = "spectrum";
		char *file_type = ".dat";
		sprintf(wfilename, "%s%d%s", base_wfilename, file_count, file_type);
		wfp = fopen(wfilename,"w"); //ファイルのオープン

		// フーリエ変換を行う
		dft(data);
		double fk;
		for (int j=0; j<DATASIZE/2; j++) {
			// printf("%d: %e\n", j, data[j]);
			fk = j / (DATASIZE * DELTA);
			fprintf(wfp,"%f %f\n", fk, data[j]) ;
		}
		/* ファイルのクローズ */
		fclose(fp);
		file_count++;
	 } // end of loop extracting files

	closedir(dir);
	return 0;	
/*
	FILE *fp;
	char *filename = argv[1];
	char readline[DATASIZE] = {'\0'};
	// ファイルのオープン
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "%sのオープンに失敗しました.\n", filename);
	}
	int i=0;
	double data[DATASIZE];
	// ファイルの終端まで文字を読み取り表示する
	// while ( fgets(readline, DATASIZE, fp) != NULL ) {
	while( fscanf(fp,"%lf",&data[i]) != EOF ){
		// printf("data[%d] = %f\n", i, data[i]);
		i++;
	}
	FILE *wfp;
	char *wfilename = "spectrum.dat";
	wfp = fopen(wfilename,"w"); //ファイルのオープン
	// フーリエ変換を行う
	dft(data);
	double fk;
	for (int j=0; j<DATASIZE/2; j++) {
		// printf("%d: %e\n", j, data[j]);
		fk = j / (DATASIZE * DELTA);
		fprintf(wfp,"%f %f\n", fk, data[j]) ;
	}
	// ファイルのクローズ
	fclose(fp);
	return 0;
*/
}
