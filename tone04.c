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
     short waveFormatType;
     short channel;
     int samplePerSec;
     int bytesPerSec;
     short blockSize;
     short bitsPerSample;
     char chunkID_data[4];
     int chunkSize_data;
     }wavehead;

//waveファイルのヘッダー作成
 void make_wave_header(wavehead *wh, int count){
     wh->chunkID_riff[0] = 'R';
     wh->chunkID_riff[1] = 'I';
     wh->chunkID_riff[2] = 'F';
     wh->chunkID_riff[3] = 'F';
     wh->chunkSize_riff = count  * 2 + 36;
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
     wh->chunkSize_data = count * 2 * wh->channel;
     }

int main(void){

 //シミュレーションのパラメータ
  double dt = (1.0 / (4.41e4)) / 50.0, dx = 0.0065; //時間、空間刻み幅
  double ts = 0.0, te = 3.0;                    //時間のスタートとエンド
  double L = 0.65;                              //弦の長さ
  int Nt = (int)((te - ts) / dt);
  //printf("%d\n", Nt);
  int Nx = (int)(L / dx);
  int i, j;
  int fin = (int)(Nx * 0.8); //弦を爪引く位置
  int p = fin;               //値を出力する点
  double newU[Nx + 1];  //振幅(両端の点をそれぞれ取るので＋1する)
  double oldU[Nx + 1];
  double newV[Nx + 1];  //速度(v[0]からv[Nx])
  double oldV[Nx + 1];
 
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
 
 //wave書き出し時に使うモノ
  double temp;
  int count = 0;

 //振幅の初期条件
  for (i = 0; i <= fin; i++){
      oldU[i] = 0.00125 * i;
      }
  for (i = fin; i <= Nx; i++){
      oldU[i] = 0.1 - 0.005 * (i - (fin));
      }

 //速度の初期条件
  newV[0] = 0.0;
  oldV[0] = 0.0;

 //境界条件（固定端）
  newU[0] = 0.0;
  oldU[0] = 0.0;
  newU[Nx] = 0.0;
  oldU[Nx] = 0.0;

 //離散化して計算,50反復ごとにdatファイル書き込み
  FILE *data;
  data = fopen("tone.dat", "wb+");
 //ここから計算 
  for (j = 1; j <= Nt; j++){
      //printf("walking\n");
      for (i = 1; i <= Nx - 1; i++){
          a1 = (oldU[i+1] - 2*oldU[i] + oldU[i-1]) / pow(dx, 2.0);
          a3 = (oldV[i+1] - 2*oldV[i] + oldV[i-1]) / pow(dx, 2.0);
          if(i == 1)
            a2 = (oldU[i+2] - 4*oldU[i+1] + 6*oldU[i] - 4*oldU[i-1] + oldU[i-1]) / pow(dx, 4.0);
          else if(i == Nx-1)
            a2 = (oldU[i+1] - 4*oldU[i+1] + 6*oldU[i] - 4*oldU[i-1] + oldU[i-2]) / pow(dx, 4.0);
          else
            a2 = (oldU[i+2] - 4*oldU[i+1] + 6*oldU[i] - 4*oldU[i-1] + oldU[i-2]) / pow(dx, 4.0);
          newV[i] = oldV[i] + k*(T*a1 - E*M*a2 - d1*oldU[i] + d3*a3);
          newU[i] = oldU[i] + dt*newV[i];
          }
      if(j%50 == 0){
          if(abs_max < fabs(newU[p])){
              abs_max = fabs(newU[p]); //ついでにabs_maxの設定
              }
          fwrite(&newU[p], sizeof(double), 1, data); //ファイル書き込み
          count++;
          }
      for (i = 0; i <= Nx; i++){
          oldU[i] = newU[i];
          oldV[i] = newV[i];
      }
  }
  //printf("swimming\n");
  fseek(data, 0, SEEK_SET); //ファイルの先頭に戻る
  //printf("%d\n", count);
 
 //waveファイルのヘッダー作成
  wavehead wh;
  make_wave_header(&wh, count);


 //ファイルに書き出す
  FILE *wave;
  wave = fopen("tone.wav", "wb");
  fwrite(&wh, sizeof(wh), 1, wave);
  for(j = 1; j <= count; j++){
      //fscanf(data, "%lf", &temp);
      fread(&temp, sizeof(double), 1, data);
      //printf("%f\n", temp);
      dat_normalized = (short)(temp / abs_max * 32767);
      fwrite(&dat_normalized, 2, 1, wave);
      //printf("running\n");
      }
  fclose(data);
  fclose(wave);

  return 0;
}
