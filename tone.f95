program tone
  implicit none

  !パラメータ設定
  real, parameter :: dt = (1/(4.41e4))/50, dx = 0.065 !時間、空間刻み幅
  real, parameter :: ts = 0.0, te = 2.0 !時間のスタートとエンド
  real, parameter :: L = 0.65 !弦の長さ
  integer, parameter :: Nt = int((te - ts)/dt)
  integer, parameter :: Nx = int(L/dx)
  integer :: i, j
  integer, parameter :: fin = int(Nx * 0.8) !弦を爪引く位置
  integer, parameter :: p = fin !値を出力する点
  real :: u(0:Nx, 0:Nt) !振幅
  real :: v(0:Nx, 0:Nt) !速度

  !弦のパラメータ
  real(8), parameter :: rho = 1140 !弦の密度
  real(8), parameter :: A = 0.5188e-6 !弦の断面積
  real(8), parameter :: T = 60.97 !弦の張力
  real(8), parameter :: E = 5.4e9 !弦のヤング率
  real(8), parameter :: M = 0.171e-12 !弦の曲げモーメント
  real(8), parameter :: d1 = 0.81e-3 !空気抵抗による振動の減衰係数
  real(8), parameter :: d3 = 0.14e-3 !粘弾性による振動の減衰係数
  real(8) :: a1, a2, a3
  real(8) :: k = dt/(rho*A)  

  !規格化とwaveファイルのパラメータ
  real :: abs_max = 0
  integer(2) :: dat_normalized
  Type wavehead
     character(4) :: chunkID_riff
     integer(4) :: chunkSize_riff
     character(4) :: formType
     character(4) :: chunkID_fmt
     integer(4) :: chunkSize_fmt
     integer(2) :: waveFormatType
     integer(2) :: channel
     integer(4) :: samplePerSec
     integer(4) :: bytesPerSec
     integer(2) :: blockSize
     integer(2) :: bitsPerSample
     character(4) :: chunkID_data
     integer(4) :: chunkSize_data
  END Type wavehead

  !waveファイルのヘッダー作成
  Type(wavehead) :: wh
  wh%chunkID_riff = "RIFF"
  wh%chunkSize_riff = (Nt + 1) * 2 + 36
  wh%formType = "WAVE"
  wh%chunkID_fmt = "fmt "
  wh%chunkSize_fmt = 16
  wh%waveFormatType = 1
  wh%channel = 1
  wh%samplePerSec = 44100
  wh%bytesPerSec = wh%samplePerSec * 2 * wh%channel
  wh%blockSize = 2 * wh%channel
  wh%bitsPerSample = 16
  wh%chunkID_data = "data"
  wh%chunkSize_data = Nt * 2 * wh%channel
  
  !振幅の初期条件 
  do i = 0, fin
     u(i, 0) = 0.00125 * i
  end do
  do i = fin, Nx
     u(i, 0) = 0.1 - 0.005*(i - (fin))
     !print *, u(i,0)
  end do

  !速度の初期条件
  do j = 0, Nt
     v(0, j) = 0
  end do  
  
  !境界条件(固定端の場合)
  do j = 0, Nt
     u(0, j) = 0.0
     u(Nx, j) = 0.0
  end do

  !print *, k
  
  !離散化して計算
  do j = 1, Nt
     do i = 1, Nx - 1
        a1 = (u(i+1, j-1) - 2*u(i, j-1) + u(i-1, j-1)) / (dx**2)
        a3 = (v(i+1, j-1) - 2*v(i, j-1) + v(i-1, j-1)) / (dx**2)
        select case(i)
        case(1)
           a2 = (u(i+2, j-1) - 4*u(i+1, j-1) + 6*u(i, j-1) - 4*u(i-1, j-1) + u(i-1, j-1)) / (dx**4)
        case(Nx - 1)
           a2 = (u(i+1, j-1) - 4*u(i+1, j-1) + 6*u(i, j-1) - 4*u(i-1, j-1) + u(i-2, j-1)) / (dx**4)
        case default
           a2 = (u(i+2, j-1) - 4*u(i+1, j-1) + 6*u(i, j-1) - 4*u(i-1, j-1) + u(i-2, j-1)) / (dx**4)
        end select        
        v(i, j) = v(i, j-1) + k * (T*a1 - E*M*a2 - d1*u(i, j-1) + d3*a3)
        u(i, j) = u(i, j-1) + dt*v(i, j)
        !print *, T*a1, E*M*a2, d3*a3
     end do  
  end do
  
  !ファイルに書き出す
  open(1, file = "tone.dat", status = "replace")
  open(2, file = "tone.wav", form = "unformatted", access = "stream")
  write(2) wh

  do j = 0, Nt
     if(mod(j,50) == 0) then
        write(1, *) j*dt, u(p, j), v(p, j)
        if(abs_max < abs(u(p, j)))then
           abs_max = abs(u(p, j))
        end if
     end if
  end do
  do j = 0, Nt
     if(mod(j,50) == 0) then
        dat_normalized = u(p, j)/abs_max*32767
        write(2) dat_normalized
     end if
  end do
  
  close(1)
  close(2)
  
end program tone
