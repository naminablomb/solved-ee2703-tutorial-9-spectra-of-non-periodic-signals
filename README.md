Download Link: https://assignmentchef.com/product/solved-ee2703-tutorial-9-spectra-of-non-periodic-signals
<br>
In the previous assignment we looked at functions that were periodic and extracted their spectra. The approach was:

<ul>

 <li>Sample the signal so that <em>f<sub>Nyquist </sub></em>is met, and so that ∆<em>f </em>is small enough. Generate the frequency axis from −<em>f<sub>max</sub></em><em>/</em>2 to +<em>f<sub>max</sub></em><em>/</em>2, taking care to drop the last term.</li>

 <li>Ensure that the signal starts at <em>t </em>= 0<sup>+ </sup>and ends at <em>t </em>= 0<sup>−</sup></li>

 <li>Use 2<em><sup>k </sup></em>samples</li>

 <li>Obtain the DFT. Rotate the samples so that they go from <em>f </em>= −<em>f<sub>max</sub></em><em>/</em>2 to <em>f </em>= +<em>f<sub>max</sub></em><em>/</em>2−∆<em>f</em>.</li>

 <li>Plot the magnitude and phase of the spectrum. Usually we would plot the magnitude in dB and the phase in degrees and the frequency axis would be logarithmic. This is to capture polynomial decay of the spectrum.</li>

</ul>

Now we want to look at non-periodic signals. Our first signal will be sin. Obtained over 0 to 2<em>π </em>with 64 samples, this function looks as follows:

<ul>

 <li>h<em>* </em>1i≡</li>

</ul>

h<em>eg1 </em>2i

<ul>

 <li>h<em>eg1 </em>2i≡ (1) from pylab import * t=linspace(-pi,pi,65);t=t[:-1] dt=t[1]-t[0];fmax=1/dt y=sin(sqrt(2)*t) y[0]=0 # the sample corresponding to -tmax should be set zeroo y=fftshift(y) # make y start with y(t=0) Y=fftshift(fft(y))/64.0 w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1] figure() subplot(2,1,1) plot(w,abs(Y),lw=2) xlim([-10,10]) ylabel(r”$|Y|$”,size=16) title(r”Spectrum of $sinleft(sqrt{2}tright)$”) grid(True) subplot(2,1,2) plot(w,angle(Y),’ro’,lw=2) xlim([-10,10]) ylabel(r”Phase of $Y$”,size=16) xlabel(r”$omega$”,size=16) grid(True) savefig(“fig10-1.png”) show()</li>

</ul>

We expected two spikes, but what we got were two peaks each with two values and a gradually decaying magnitude. The phase is correct though – but take a look at the code. There is one line there:

y[0]=0 # the sample corresponding to -tmax should be set zero

What is this for? An antisymmetric function has a purely imaginary fourier transform. But what about an antisymmetric set of samples? Suppose

<em>y</em>[0] =            0            sin(0)

<em>N</em>

<em>y</em>[<em>i</em>] =      −<em>y</em>[<em>N </em>−<em>i</em>]       <em>i </em>= 1<em>,</em>2<em>…    </em> −1

2

<em>N</em>

<em>y</em><em> t</em><em>N</em>

The DFT of this sequence will give us

<em>Y</em>[<em>k</em>]     =

<em>n</em>=0                                <em>N</em>

<em>N</em>

<em>N</em>

=               <em>y</em>[<em>n</em>]    exp    <em>kn     </em>−exp     − <em>kn       </em>+<em>y</em>[  ]exp(<em>πkj</em>)

=1                                       <em>N                               N                      </em>2

<em>N</em>

<em><sub>k       </sub>N</em>

=                                                       +(−1) <em>y</em>[  ]

2

But this is no longer pure imaginary! Since we need that property, we set <em>y</em>[<em>N</em><em>/</em>2] to zero. That gives us a purely imaginary phase. But what to do about the magnitude? And what went wrong?

To understand what went wrong, let us plot the time function over several time periods.

<ul>

 <li>h<em>eg2 </em>4i≡ from pylab import * t1=linspace(-pi,pi,65);t1=t1[:-1] t2=linspace(-3*pi,-pi,65);t2=t2[:-1] t3=linspace(pi,3*pi,65);t3=t3[:-1]</li>

</ul>

# y=sin(sqrt(2)*t) figure(2) plot(t1,sin(sqrt(2)*t1),’b’,lw=2) plot(t2,sin(sqrt(2)*t2),’r’,lw=2) plot(t3,sin(sqrt(2)*t3),’r’,lw=2) ylabel(r”$y$”,size=16) xlabel(r”$t$”,size=16) title(r”$sinleft(sqrt{2}tright)$”) grid(True) savefig(“fig10-2.png”) show()

The blue line connects the points whose DFT we took. The red lines show the continuation of the function. Quite clearly, even though sin is a periodic function, the portion between −<em>π </em>and <em>π </em>is not the part that can be replicated to generate the function. So which function is the DFT trying to fourier analyse? For that we have to replicate just the blue points. And that is shown below:

<ul>

 <li>h<em>eg3 </em>5i≡ from pylab import * t1=linspace(-pi,pi,65);t1=t1[:-1] t2=linspace(-3*pi,-pi,65);t2=t2[:-1] t3=linspace(pi,3*pi,65);t3=t3[:-1] y=sin(sqrt(2)*t1) figure(3) plot(t1,y,’bo’,lw=2) plot(t2,y,’ro’,lw=2) plot(t3,y,’ro’,lw=2) ylabel(r”$y$”,size=16) xlabel(r”$t$”,size=16)</li>

</ul>

title(r”$sinleft(sqrt{2}tright)$ with $t$ wrapping every $2pi$ “) grid(True) savefig(“fig10-3.png”) show()

<h1>Gibbs Phenomenon</h1>

This is something you know very well. The Fourier transform of the box function

<em>f</em>

is given by

2sin(<em>ωt</em><sub>0</sub>)

<em>F</em>(<em>ω</em>) = <em>ω</em>

The spectrum of the box function decays very slowly, as 2<em>/ω</em>.

Now our function is an odd function with a big jump. So let us consider the periodic ramp:

<em>f</em>(<em>t</em>) = <em>t</em><em>, </em>−<em>π </em><em>&lt; t </em><em>&lt; π</em>

Then the fourier series of this ramp is

<em>f</em>

Again the coefficients decay very slowly.

The DFT is just like the fourier series, except that both time and frequency are samples. So, if the time samples are like a ramp, the frequency samples will decay as 1<em>/ω</em>. Let us verify this for the ramp itself

7 h<em>eg4 </em>7i≡ from pylab import * t=linspace(-pi,pi,65);t=t[:-1] dt=t[1]-t[0];fmax=1/dt y=t y[0]=0 # the sample corresponding to -tmax should be set zeroo y=fftshift(y) # make y start with y(t=0) Y=fftshift(fft(y))/64.0 w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1] figure() semilogx(abs(w),20*log10(abs(Y)),lw=2) xlim([1,10]) ylim([-20,0]) xticks([1,2,5,10],[“1″,”2″,”5″,”10″],size=16) ylabel(r”$|Y|$ (dB)”,size=16) title(r”Spectrum of a digital ramp”) xlabel(r”$omega$”,size=16) grid(True) savefig(“fig10-4.png”) show()

Clearly the spectrum decays as 20 dB per decade, which corresponds to 1<em>ω</em>. The big jumps at <em>nπ </em>force this slowly decaying spectrum, which is why we don’t see the expected spikes for the spectrum of sin.

<h1>Windowing</h1>

So what do we do? Well the spikes happen at the end of the periodic interval. So we damp the function near there, i.e., we multiply our function sequence <em>f</em>[<em>n</em>] by a “window” sequence <em>w</em>[<em>n</em>]: <em>g</em>(<em>n</em>) = <em>f</em>(<em>n</em>)<em>w</em>(<em>n</em>)

The new spectrum is got by convolving the two fourier transforms:

<h2><em>N</em>−1</h2>

<em>G<sub>k </sub></em>= ∑ <em>F</em><em>nW<sub>k</sub></em>−<em><sub>n </sub></em><em>n</em>=0

Suppose <em>f<sub>n </sub></em>is a sinusoid. Then <em>F<sub>k </sub></em>has two spikes. But the two spikes are now smeared out by <em>W<sub>k</sub></em>. So we expect to get broader peaks. But what this also does is to suppress the jump at the edge of the window. The window we will use is called the Hamming window:

<em>w</em>

Let us look at our time sequence for sin now …

8              h<em>eg5 </em>8i≡ from pylab import * t1=linspace(-pi,pi,65);t1=t1[:-1] t2=linspace(-3*pi,-pi,65);t2=t2[:-1] t3=linspace(pi,3*pi,65);t3=t3[:-1] n=arange(64) wnd=fftshift(0.54+0.46*cos(2*pi*n/63)) y=sin(sqrt(2)*t1)*wnd figure(3) plot(t1,y,’bo’,lw=2) plot(t2,y,’ro’,lw=2) plot(t3,y,’ro’,lw=2) ylabel(r”$y$”,size=16) xlabel(r”$t$”,size=16) title(r”$sinleft(sqrt{2}tright)times w(t)$ with $t$ wrapping every $2pi$ “) grid(True) savefig(“fig10-5.png”) show()

The jump is still there, but it is much reduced. There is a little bit of magic in keeping some of the jump – it gives us an extra 10 db of suppression.

Now let us take the DFT of this sequence and see what we get:

10 h<em>eg6 </em>10i≡ from pylab import * t=linspace(-pi,pi,65);t=t[:-1] dt=t[1]-t[0];fmax=1/dt n=arange(64) wnd=fftshift(0.54+0.46*cos(2*pi*n/63)) y=sin(sqrt(2)*t)*wnd y[0]=0 # the sample corresponding to -tmax should be set zeroo y=fftshift(y) # make y start with y(t=0) Y=fftshift(fft(y))/64.0

w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1] figure() subplot(2,1,1) plot(w,abs(Y),lw=2) xlim([-8,8]) ylabel(r”$|Y|$”,size=16) title(r”Spectrum of $sinleft(sqrt{2}tright)times w(t)$”) grid(True) subplot(2,1,2) plot(w,angle(Y),’ro’,lw=2) xlim([-8,8])

ylabel(r”Phase of $Y$”,size=16) xlabel(r”$omega$”,size=16) grid(True) savefig(“fig10-6.png”)

show()

Compare to our first plot and you can see that the magnitude is greatly improved. We√ still have a peak that is two samples wide. But that is because 2 lies between 1 and 2, which are the two fourier components available. If we use four times the number of points we should get better results.

12 h<em>eg7 </em>12i≡ from pylab import * t=linspace(-4*pi,4*pi,257);t=t[:-1] dt=t[1]-t[0];fmax=1/dt n=arange(256) wnd=fftshift(0.54+0.46*cos(2*pi*n/256)) y=sin(sqrt(2)*t) # y=sin(1.25*t) y=y*wnd

y[0]=0 # the sample corresponding to -tmax should be set zeroo y=fftshift(y) # make y start with y(t=0) Y=fftshift(fft(y))/256.0 w=linspace(-pi*fmax,pi*fmax,257);w=w[:-1] figure() subplot(2,1,1) plot(w,abs(Y),’b’,w,abs(Y),’bo’,lw=2) xlim([-4,4]) ylabel(r”$|Y|$”,size=16) title(r”Spectrum of $sinleft(sqrt{2}tright)times w(t)$”) grid(True) subplot(2,1,2) plot(w,angle(Y),’ro’,lw=2) xlim([-4,4]) ylabel(r”Phase of $Y$”,size=16)

xlabel(r”$omega$”,size=16) grid(True) savefig(“fig10-7.png”) show()

Is it better? Well it is quite a bit better since we are now zoomed in and see a lot more detail. But why is it not just a single peak? The reason for that is <em>w</em>(<em>t</em>). Multiplication in time is convolution in frequency and vice versa. So by multiplying with <em>w</em>(<em>t</em>), we got rid of the 1<em>/f </em>decay. But the delta function is now replaced by the shape of the DFT of <em>w</em>[<em>n</em>]. That gives us a factor of two broadening over the peak when there is no window, which is why we still see a peak whose width is two samples.√

Note that it is <em>not </em>because 2 is between 1<em>.</em>25 and 1<em>.</em>5. To verify, there is an alternate function in the above code, namely sin(1<em>.</em>25<em>t</em>). But this gives a broad peak as well. That is because of <em>w</em>[<em>n</em>].

<h1>The Assignment</h1>

<ol>

 <li>Work through the example codes and understand them.</li>

 <li>Consider the function cos<sup>3</sup>(<em>ω</em><sub>0</sub><em>t</em>). Obtain its spectrum for <em>ω</em><sub>0 </sub>= 0<em>.</em>86 with and without a Hamming window.</li>

 <li>Write a program that will take a 128 element vector known to contain cos(<em>ω</em><sub>0</sub><em>t </em>+<em>δ</em>) for arbitrary <em>δ </em>and 0<em>.</em>5 <em>&lt; ω</em><sub>0 </sub><em>&lt; </em>1<em>.</em> The values of <em>t </em>go from −<em>π </em>to <em>π</em>. You have to extract the digital spectrum of the signal, find the two peaks at ±<em>ω</em><sub>0</sub>, and estimate <em>ω</em><sub>0 </sub>and <em>δ</em>.</li>

 <li>Suppose the data in Q3 has added “white gaussian noise”. This can be generated byrandn() in python. The extent of this noise is 0<em>.</em>1 in amplitude (i.e., 0.1*randn(N), where N is the number of samples). Repeat the problem and find the <em>ω</em><sub>0 </sub>and <em>δ</em></li>

 <li>Plot the DFT of the function</li>

</ol>

for <em>t </em>going from −<em>π </em>to <em>π </em>in 1024 steps. This is known as a “chirped” signal, and its frequency continuously changes from 16 to 32 radians per second. This also means that the period is 64 samples near −<em>π </em>and is 32 samples near +<em>π</em>.

<ol start="6">

 <li>For the same chirped signal, break the 1024 vector into pieces that are 64 sampleswide. Extract the DFT of each and store as a column in a 2D array. Then plot the array as a surface plot to show how the frequency of the signal varies with</li>

</ol>

This is new. So far we worked either in time or in frequency. But this is a “timefrequency” plot, where we get localized DFTs and show how the spectrum evolves in time.