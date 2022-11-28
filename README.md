Levenberg-Marquardt nonlinear least squares algorithms for C# (.NET users) for x64

Base calculation code from http://users.ics.forth.gr/~lourakis/levmar/index.html

Code base on https://github.com/AvengerDr/LevmarSharp but this version support only x64 (not both)

Levmar side is same signature;

<code> LevmarDif, LevmarDer </code>

Added convolution&deconvolution based on fft transform;

FFT library's from http://www.fftw.org/download.html

  <code> fft, stepresponse_conv, stepresponse_deconv, conv ,deconv  </code>

For usage; have to use both dll (fft library not referancable, but have to in same directory)
levmarSharp_x64.dll
libfftw3-3.dll
