Levenberg-Marquardt nonlinear least squares algorithms for C# (.NET users) for x64

Base calculation code from http://users.ics.forth.gr/~lourakis/levmar/index.html

Code base on https://github.com/AvengerDr/LevmarSharp but this version support only x64 (not both)

Levmar side is same signature;

<code>
int LevmarDif(LevmarFunc^ f, array<double>^% p, array<double>^ x, int m, int n, int itmax, array<double>^ opts, array<double>^% info, IntPtr data);
  
int LevmarDif(LevmarFunc^ f, array<double>^% p, array<double>^ x, int m, int n, int itmax, array<double>^ opts, array<double>^% info);

int LevmarDer(LevmarFunc^ f, LevmarJacf^ j, array<double>^% p, array<double>^ x, int m, int n, int itmax, array<double>^ opts, array<double>^% info);
  
int LevmarDer(LevmarFunc^ f, LevmarJacf^ j, array<double>^% p, array<double>^ x, int m, int n, int itmax, array<double>^ opts, array<double>^% info, IntPtr data);
  </code>

Added convolution&deconvolution based on fft transform;

FFT library's from http://www.fftw.org/download.html

  <code>
int fft(array<double>^% out, array<double>^ in, int n);

int stepresponse_conv_integral(array<double>^% g, array<double>^ x, array<double>^ h, int n, int m );
    
int stepresponse_conv(array<double>^% g, array<double>^ x, array<double>^ h, int n, int m);
    
int stepresponse_deconv(array<double>^% g, array<double>^ x, array<double>^ h, int n, int m);
    

int conv(array<double>^% g, array<double>^ x, array<double>^ h, int n, int m);
    
int deconv(array<double>^% g, array<double>^ x, array<double>^ h, int n, int m);
  </code>

For usage; have to use both dll (fft library not referancable, but have to in same directory)
levmarSharp_x64.dll
libfftw3-3.dll
