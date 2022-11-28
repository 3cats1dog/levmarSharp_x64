<p>Levenberg-Marquardt nonlinear least squares algorithms for C# (.NET users) for x64</p>
<p>Base calculation code from http://users.ics.forth.gr/~lourakis/levmar/ without LAPACK</p>
<br />
<p>Dll wrapper based on https://github.com/AvengerDr/LevmarSharp but this version support only x64 (not both)</p>
<p><i>* The C# wrapper was implemented by Adalberto L. Simeone (http://www.adalsimeone.me) and is released under the terms of the GPLv3 license.) </i>

<p>Levmar side is same signature with <a href="https://github.com/AvengerDr/LevmarSharp">LevmarSharp</a>;</p>

<code> LevmarDif, LevmarDer </code>

<p>Added convolution and deconvolution based on fft transform;</p>
<code> fft, conv ,deconv, stepresponse_conv, stepresponse_deconv </code>
<br/>
<br/>
FFT library's from http://www.fftw.org/download.html
<br/>
<p><b>For usage;</b> have to use both dll (fft library not referenceable, but have to in same directory)</p>
<ul>
  <li>levmarSharp_x64.dll</li>
  <li>libfftw3-3.dll</li>
</ul>

