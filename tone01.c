#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//waveファイルのヘッダーのひな型
typedef struct{
 char chunkID_riff[4];
 int chunkSize_riff;
 char formType[4];
 char chunkID_fmt[4];
 int chunkSize_fmt;
 short int waveFormatType;
 short int channel;
 int samplePerSec;
 int bytesPerSec;
 short int blockSize;
 short int bitsPerSample;
 char chunkID_data[4];
 int chunkSize_data;
 }wavehead;

//waveファイルのヘッダー作成
 void make_wave_header(wavehead *wh, int Nt){
     wh->chunkID_riff[0] = 'R';
     wh->chunkID_riff[1] = 'I';
     wh->chunkID_riff[2] = 'F';
     wh->chunkID_riff[3] = 'F';
     wh->chunkSize_riff = (Nt + 1) * 2 + 36;
     wh->formType[0] = 'W';
     wh->formType[1] = 'A';
     wh->formType[2] = 'V';
     wh->formType[3] = 'E';
     wh->chunkID_fmt[0]='f';
     wh->chunkID_fmt[1]='m';
     wh->chunkID_fmt[2]='t';
     wh->chunkID_fmt[3]=' ';
     wh->chunkSize_fmt = 16;
     wh->waveFormatType = 1;
     wh->channel = 1;
     wh->samplePerSec = 44100;
     wh->bytesPerSec = wh->samplePerSec * 2 * wh->channel;
     wh->blockSize = 2 * wh->channel;
     wh->bitsPerSample = 16;
     wh->chunkID_data[0] = 'd';
     wh->chunkID_data[1] = 'a';
     wh->chunkID_data[2] = 't';
     wh->chunkID_data[3] = 'a';
     wh->chunkSize_data = Nt * 2 * wh->channel;
 }

int main(void){
    printf("running\n");
    //シミュレーションのパラメータ
    double dt = (1.0 / (4.41e4)) / 50.0, dx = 0.065; //時間、空間刻み幅
    printf("%e\n",dt);
    double ts = 0.0, te = 2.0;                    //時間のスタートとエンド
    double L = 0.65;                              //弦の長さ
    int Nt = (int)((te - ts) / dt);
    printf("%d\n", Nt);
    int Nx = (int)(L / dx);
    int i, j;
    int fin = (int)(Nx * 0.8); //弦を爪引く位置
    int p = fin;               //値を出力する点
    printf("running1\n");
    double u[Nx + 1][Nt + 1];  //振幅(両端の点をそれぞれ取るので＋1する)
    printf("running3\n");
    fflush(stdout);
    double v[Nx + 1][Nt + 1];  //速度(v[0]からv[Nt])
    printf("runnning2\n");
    //弦のパラメータ（ここをいじると音色が変わります）
    double rho = 1140;    //弦の密度
    double A = 0.5188e-6; //弦の断面積
    double T = 60.97;     //弦の張力
    double E = 5.4e9;     //弦のヤング率
    double M = 0.171e-12; //弦の曲げモーメント
    double d1 = 0.81e-3;  //空気抵抗による振動の減衰係数
    double d3 = 0.14e-3;  //粘弾性による振動の減衰係数
    double a1, a2, a3;
    double k = dt / (rho * A);

    //規格化のパラメータ
    double abs_max = 0;
    short int dat_normalized;

    //waveファイルのヘッダー作成
    wavehead wh;
    make_wave_header(&wh, Nt);

    //振幅の初期条件
    for (i = 0; i <= fin; i++)
    {
        u[i][0] = 0.00125 * i;
  }
  for (i = fin; i <= Nx; i++){
      u[i][0] = 0.1 - 0.005 * (i - (fin));
  }

 //速度の初期条件
  for (j = 0; j <= Nt; j++){
      v[0][j] = 0.0;
  }

 //境界条件（固定端）
  for (j = 0; j <= Nt; j++){
      u[0][j] = 0.0;
      u[Nx][j] = 0.0;
  }

 //離散化して計算
  for (j = 1; j <= Nt; j++){
      for (i = 1; i <= Nx - 1; i++){
          a1 = (u[i+1][j-1] - 2*u[i][j-1] + u[i-1][j-1]) / pow(dx, 2.0);
          a3 = (v[i+1][j-1] - 2*v[i][j-1] + v[i-1][j-1]) / pow(dx, 2.0);
          if(i == 1)
            a2 = (u[i+2][j-1] - 4*u[i+1][j-1] + 6*u[i][j-1] - 4*u[i-1][j-1] + u[i-1][j-1]) / pow(dx, 4.0);
          else if(i == Nx - 1)
            a2 = (u[i+1][j-1] - 4*u[i+1][j-1] + 6*u[i][j-1] - 4*u[i-1][j-1] + u[i-2][j-1]) / pow(dx, 4.0);
          else
            a2 = (u[i+2][j-1] - 4*u[i+1][j-1] + 6*u[i][j-1] - 4*u[i-1][j-1] + u[i-2][j-1]) / pow(dx, 4.0);
          v[i][j] = v[i][j-1] + k*(T*a1 - E*M*a2 - d1*u[i][j-1] + d3*a3);
          u[i][j] = u[i][j - 1] + dt*v[i][j];
      }
      //ついでなのでabs_maxもここで設定する 
      if(j%50 == 0){
          if(abs_max < fabs(u[p][j])){
              abs_max = fabs(u[p][j]);
          }
      }
  }

 //ファイルに書き出す
  FILE *file;
  file = fopen("c.wav", "wb");
  fwrite(&wh, sizeof(wh), 1, file);
  for (j = 0; j <= Nt; j++){
      if (j % 50 == 0){
          dat_normalized = (short)(u[p][j] / abs_max * 32767);
          fwrite(&dat_normalized, 2, 1, file);
      }
  }
  fclose(file);
  

  return 0;
}
