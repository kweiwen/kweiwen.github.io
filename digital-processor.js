function comparator(input_data_1, input_data_2) { return (input_data_1 > input_data_2) ? 1 : (input_data_1 > input_data_2) ? -1 : 0; }

function fft(re, im) 
{
  var N = re.length;
  for (var i = 0; i < N; i++) 
  {
    for(var j = 0, h = i, k = N; k >>= 1; h >>= 1)
      j = (j << 1) | (h & 1);
    if (j > i) 
    {
      re[j] = [re[i], re[i] = re[j]][0]
      im[j] = [im[i], im[i] = im[j]][0]
    }
  }
  for(var hN = 1; hN * 2 <= N; hN *= 2)
  {
    for (var i = 0; i < N; i += hN * 2)
    {
      for (var j = i; j < i + hN; j++) 
      {
        var cos = Math.cos(Math.PI * (j - i) / hN),
        sin = Math.sin(Math.PI * (j - i) / hN)
        var tre =  re[j+hN] * cos + im[j+hN] * sin,
        tim = -re[j+hN] * sin + im[j+hN] * cos;
        re[j + hN] = re[j] - tre; im[j + hN] = im[j] - tim;
        re[j] += tre; im[j] += tim;
      }
    }
  }
}

function ifft(re, im)
{
  fft(im, re);
  for(let i = 0, N = re.length; i < N; i++)
  {
    im[i] = im[i] / N;
    re[i] = re[i] / N;
  }
}

function polar2complex(amplitude, angle)
{
  let real = amplitude * Math.cos(angle);
  let imag = amplitude * Math.sin(angle);
  return [real, imag];
}

function roll(array, k)
{
  let step = k % array.length;
  for(let i = 0 ; i < step ; i++)
  {
      let value = array.pop();
      array.unshift(value);
  }
}

function complex_roll(array, k)
{
  let size = array.length;
  let half = size / 2;
  
  let dc = array[0];
  let arr1 = array.slice(1, half);
  let nyquist = array[half];
  let arr2 = array.slice(half + 1, size);
  
  roll(arr1, k);
  roll(arr2, half - k - 1);
  
  let result = [];
  result = result.concat([dc], arr1, [nyquist], arr2);
  return result;
}

function sigmoid(input_data, gain)
{
  input_data = input_data * gain;
  return input_data / Math.sqrt(input_data * input_data + 1);
}

function gaussian()
{
  let R1 = Math.random();
  let R2 = Math.random();
  return  Math.sqrt(Math.log(R1) * -2) * Math.cos(2 * Math.PI * R2);
}

function clip(num, min, max) 
{
  return num <= min 
         ? min 
         : num >= max 
         ? max 
         : num
}

function vca(voltage) 
{ 
  return Math.pow(10, -1.5 * voltage); 
}

class lorenz_attractor
{
  constructor(options)
  {
    this.sigma = 10;
    this.rho = 13.927;
    this.beta = 8/3;
    this.time = 0.01;
    
    this.x = Math.random();
    this.y = Math.random();
    this.z = Math.random();

    this.flag = 0;
  }
  
  lorenz_attr(x, y, z)
  {
    let x_dot = this.sigma * (y - x);
    let y_dot = this.rho * x - y - x * z;
    let z_dot = x * y - this.beta * z;
    return [x_dot, y_dot, z_dot];
  }
  
  process()
  {
    let result = this.lorenz_attr(this.x, this.y, this.z);
    let x_dot = result[0];
    let y_dot = result[1];
    let z_dot = result[2];
    this.x = this.x + (x_dot * this.time);
    this.y = this.y + (y_dot * this.time);
    this.z = this.z + (z_dot * this.time);
  }

  sample(clock)
  {
    if(clock == 1 & this.flag == 0)
    {
      this.process();
      this.flag = 1;
    }
    else if(clock == 0 & this.flag == 1)
    {
      this.flag = 0;
    }
  }
}

class shift_register
{
  constructor(options)
  {
    this.flag = 0;
    this.buffer = new Array(8).fill(-1);
    this.index = 0;
    this.mask = 7;
  }

  sample(clock, input_data)
  {
    if(clock > 0 & this.flag == 0)
    {
      this.buffer[this.index] = input_data
      this.index = (this.index + 1) & this.mask
      this.flag = 1
    }
    else if(clock <= 0 & this.flag == 1)
    {
      this.flag = 0
    }
  }
  
  xor(p1, q8)
  {
    // P2(CLK)
    // Q8       L   H   L   H
    // P1       L   L   H   H
    // Output  -1   1   1  -1

    if(q8 < 0)
    {
      return p1;
    }
    else if(q8 > 0)
    {
      return -p1;
    }
    else
    {
      return 0;
    }
  }
}

class svf
{
  constructor(options, cutoffFrequency, resonance)
  {
    this.sampleRate = sampleRate;
    this._invSampleRate = 1 / this.sampleRate;
    this._cutoffFrequency = cutoffFrequency;
    this._resonance = resonance;
    this.update();
    this.reset();
  }

  update()
  {
    this.g  = (Math.PI * this._cutoffFrequency * this._invSampleRate);
    this.R2 = (1.0 / this._resonance);
    this.h  = (1.0 / (1.0 + this.R2 * this.g + this.g * this.g));
  }

  reset()
  {
    this.s1 = 0;
    this.s2 = 0;
    
    this.HP = 0;
    this.BP = 0;
    this.LP = 0;
  }

  set cutoffFrequency(input_data)
  {
    this._cutoffFrequency = input_data;
    this.update();
  }

  get cutoffFrequency()
  {
    return this._cutoffFrequency;
  }

  set resonance(input_data)
  {
    this._resonance = input_data;
    this.update();
  }

  get resonance()
  {
    return this._resonance;
  }

  process(sample)
  {
    this.HP = (sample - this.s1 * this.R2 - this.s1 * this.g - this.s2) * this.h;
    
    this.BP = this.HP * this.g + this.s1;
    this.s1 = this.HP * this.g + this.BP;
    
    this.LP = this.BP * this.g + this.s2;
    this.s2 = this.BP * this.g + this.LP;
    
    return {"HP": this.HP, "BP": this.BP, "LP": this.LP}
  }
}

class functionGenerator
{
  constructor(options)
  {
    this.sampleRate = sampleRate;
    this._invSampleRate = 1 / this.sampleRate;
    this._invPi = 1 / Math.PI;
		this._degree = 0;
		this._cycle = 0;
    this._pitch = 0;
  }

  flushSamples()
  {
    this._degree = 0;
		this._cycle = 0;
  }

  Generate()
  {
    this._cycle = this._cycle + this._pitch * this._invSampleRate;
    this._cycle = this._cycle - Math.floor(this._cycle);
    this._degree = this._cycle * Math.PI * 2;

    let triangle = 2 * this._invPi * Math.asin(Math.sin(this._degree))
    let square = (this._degree < 0.5) ? 1 : (this._degree > 0.5) ? -1 : 0;

    return {"square": square, "triangle": triangle};
  }

  GenerateSine()
  {
    this._cycle = this._cycle + this._pitch * this._invSampleRate;
    this._cycle = this._cycle - Math.floor(this._cycle);
    this._degree = this._cycle * Math.PI * 2;

    return Math.sin(this._degree);
  }
}

class circularBuffer
{
  constructor(options)
  {
    this._writeIndex = 0;
    this._bufferLength = 0;
    this._wrapMask = 0;
    this._bufferLength = 0;
    this._buffer = null;
  }

  createBuffer(input)
  {
    this._writeIndex = 0;
    this._bufferLength = parseInt(Math.pow(2, Math.ceil(Math.log10(input) / Math.log10(2))));
    this._wrapMask = this._bufferLength - 1;
    this._buffer = new Array(this._bufferLength).fill(0);
  }

  flushBuffer()
  {
    for (let index = 0; index < this._bufferLength; index++)
    {
      this._buffer[index] = 0;
    } 
  }

  writeBuffer(input)
  {
    this._buffer[this._writeIndex++] = input;
    this._writeIndex = this._writeIndex & this._wrapMask;
  }

  readBuffer(delayInSamples)
  {
    let index = this._writeIndex - delayInSamples;
    index = index & this._wrapMask;
    return this._buffer[index];
  }

  readBufferReverse(delayInSamples)
  {
    let index = this._writeIndex + delayInSamples;
    index = index & this._wrapMask;
    return this._buffer[index];
  }

  readBufferLinear(delayInSamples)
  {
    let y1 = this.readBuffer(parseInt(delayInSamples));
    let y2 = this.readBuffer(parseInt(delayInSamples) + 1);
    let fraction = delayInSamples - parseInt(delayInSamples);
    if (fraction >= 1.0) return y2;
    return fraction * y2 + (1 - fraction) * y1;
  }

  readBufferLinearReverse(delayInSamples)
  {
    let y1 = this.readBufferReverse(parseInt(delayInSamples));
    let y2 = this.readBufferReverse(parseInt(delayInSamples) + 1);
    let fraction = delayInSamples - parseInt(delayInSamples);
    if (fraction >= 1.0) return y2;
    return fraction * y2 + (1 - fraction) * y1;
  }
}

class coefficient
{
    constructor(options)
    {
      this._filterType = options.filterType; 
      this._frequencyCut = options.frequencyCut;
      this._sampleRate = options.sampleRate;
      this._res = options.res;
      this._slope = options.slope;
      this._magnitude = options.magnitude;
      
      this.EULER = Math.exp(1);

      this.a0 = 0;
      this.a1 = 0;
      this.a2 = 0;
      this.b0 = 0;
      this.b1 = 0;
      this.b2 = 0;

      this.generate();
    }

    get omega()
    {
      return (2 * Math.PI * this._frequencyCut / this._sampleRate);
    }
      
    get gain()
    {
      return Math.pow(10, (this._magnitude / 20));
    }
    
    get sineOmega()
    {
      return Math.sin(this.omega);
    }

    get cosineOmega()
    {
      return Math.cos(this.omega);
    }
    
    set filterType(input_data)
    {
      this._filterType = input_data;
      this.generate();
    }

    get filterType()
    {
      return this._filterType;
    }

    set frequencyCut(input_data)
    {
      this._frequencyCut = input_data;
      this.generate();
    }

    get frequencyCut()
    {
      return this._frequencyCut;
    }

    set res(input_data)
    {
        this._res = input_data;
        this.generate();
    }
    
    get res()
    {
        return this._res;
    }
    
    generate()
    {
      let filter_type = this._filterType;
      
      if(filter_type == 'flat')
      {
        this.b0 = 1.0;
        this.b1 = 0.0;
        this.b2 = 0.0;
        
        this.a1 = 0.0;
        this.a2 = 0.0;
      }
      else if(filter_type == 'low-pass1')
      {
        this.a1 = Math.pow(this.EULER, -this.omega);
        
        this.b0 = (1.0 - this.a1) * this.gain;
        this.b1 = 0.0;
        this.b2 = 0.0;
        
        this.a1 = -this.a1;
        this.a2 = 0.0; 
      }
      else if(filter_type == 'high-pass1')
      {
        this.a1 = Math.pow(this.EULER, -this.omega);
        
        this.b0 = (1.0 + this.a1) * this.gain / 2;
        this.b1 = -(1.0 + this.a1) * this.gain / 2;
        this.b2 = 0.0;
        
        this.a1 = -this.a1;
        this.a2 = 0.0;
      }
      else if(filter_type == 'all-pass1')
      {
        this.a1 = Math.pow(this.EULER, -this.omega);
        
        this.b0 = -this.a1 * this.gain;
        this.b1 = this.gain;
        this.b2 = 0.0;
        
        this.a1 = -this.a1;
        this.a2 = 0.0;
      }
      else if(filter_type == 'low-pass2')
      {
        let alpha = this.sineOmega / (2 * this.res);

        this.a0 = 1 + alpha;
        this.a1 = (-2 * this.cosineOmega) / this.a0;
        this.a2 = (1 - alpha) / this.a0;
        
        this.b0 = (1 - this.cosineOmega) * this.gain / 2 / this.a0;
        this.b1 = (1 - this.cosineOmega) * this.gain / this.a0;
        this.b2 = (1 - this.cosineOmega) * this.gain / 2 / this.a0;
      }
      else if(filter_type == 'high-pass2')
      {
        let alpha = this.sineOmega / (2 * this.res);
        
        this.a0 = 1 + alpha;
        this.a1 = (-2 * this.cosineOmega) / this.a0;
        this.a2 = (1 - alpha) / this.a0;
        
        this.b0 = (1 + this.cosineOmega) * this.gain / 2 / this.a0;
        this.b1 = -(1 + this.cosineOmega) * this.gain / this.a0;
        this.b2 = (1 + this.cosineOmega) * this.gain / 2 / this.a0;
      }
      else if(filter_type == 'all-pass2')
      {
        let alpha = this.sineOmega / (2 * this.res);
        
        this.a0 = 1 + alpha;
        this.a1 = (-2 * this.cosineOmega) / this.a0;
        this.a2 = (1 - alpha) / this.a0;
        
        this.b0 = (1 - alpha) * this.gain / this.a0;
        this.b1 = (-2 * this.cosineOmega) * this.gain / this.a0;
        this.b2 = (1 + alpha) * this.gain / this.a0;
      }
      else if(filter_type == 'parametric')
      { 
        let A = Math.pow(10.0, (this.magnitude / 40.0))
        let alpha = this.sineOmega / (2 * A * this.res);
        
        this.a0 = 1 + alpha / A;
        this.a1 = (-2 * this.cosineOmega) / this.a0;
        this.a2 = (1 - alpha / A) / this.a0;
        
        this.b0 = (1 + alpha * A) / this.a0;
        this.b1 = -(2 * this.cosineOmega) / this.a0;
        this.b2 = (1 - alpha * A) / this.a0;
      }
      else if(filter_type == 'peaking')
      {
        let A = Math.pow(10.0, (this.magnitude / 40.0));
        let alpha = this.sineOmega / (2 * this.res);
        
        this.a0 = 1 + alpha / A;
        this.a1 = (-2 * this.cosineOmega) / this.a0;
        this.a2 = (1 - alpha / A) / this.a0;
        
        this.b0 = (1 + alpha * A) / this.a0;
        this.b1 = -(2 * this.cosineOmega) / this.a0;
        this.b2 = (1 - alpha * A) / this.a0;
      }
      else if(filter_type == 'band-pass')
      {
        let alpha = this.sineOmega / (2 * this.res);
        
        this.a0 = 1 + alpha;
        this.a1 = (-2.0 * this.cosineOmega) / this.a0;
        this.a2 = (1 - alpha) / this.a0;
        
        this.b0 = (alpha * this.gain) / this.a0;
        this.b1 = 0.0;
        this.b2 = -(alpha * this.gain) / this.a0;
      }
      else if(filter_type == 'notch')
      {
        let alpha = this.sineOmega / (2 * this.res);
        
        this.a0 = 1 + alpha;
        this.a1 = (-2.0 * this.cosineOmega) / this.a0;
        this.a2 = (1.0 - alpha) / this.a0;
        
        this.b0 = this.gain / this.a0;
        this.b1 = (-2.0 * this.cosineOmega * this.gain) / this.a0;
        this.b2 = this.gain / this.a0;
      }
      else if(filter_type == 'low-shelf')
      {
        if (this.slope >= 2.0)
          this.slope = 2.0;
        elif (this.slope <= 0.1)
          this.slope = 0.1;
        
        let A = Math.pow(10.0, (this.magnitude / 40.0));
        let alpha = this.sineOmega / 2 * Math.sqrt((A + 1 / A) * (1 / this.slope - 1) + 2);
        
        this.a0 = (A + 1) + (A - 1) * this.cosineOmega + 2 * Math.sqrt(A) * alpha;
        this.a1 = (-2 * ((A - 1) + (A + 1) * this.cosineOmega)) / this.a0;
        this.a2 = ((A + 1) + (A - 1) * this.cosineOmega - 2 * Math.sqrt(A) * alpha) / this.a0;
        
        this.b0 = (A * ((A + 1) - (A - 1) * this.cosineOmega + 2 * Math.sqrt(A) * alpha)) / this.a0;
        this.b1 = (2 * A * ((A - 1) - (A + 1) * this.cosineOmega)) / this.a0;
        this.b2 = (A * ((A + 1) - (A - 1) * this.cosineOmega - 2 * Math.sqrt(A) * alpha)) / this.a0;
      }
      else if(filter_type == 'high-shelf')
      {

        if (this.slope >= 2.0)
          this.slope = 2.0;
        elif (this.slope <= 0.1)
          this.slope = 0.1;
        
        let A = Math.pow(10.0, (this.magnitude / 40.0));
        let alpha = this.sineOmega / 2 * Math.sqrt((A + 1 / A) * (1 / this.slope - 1) + 2);
        
        this.a0 = (A + 1) - (A - 1) * this.cosineOmega + 2 * Math.sqrt(A) * alpha;
        this.a1 = (2 * ((A - 1) - (A + 1) * this.cosineOmega)) / this.a0;
        this.a2 = ((A + 1) - (A - 1) * this.cosineOmega - 2 * Math.sqrt(A) * alpha) / this.a0;
        
        this.b0 = (A * ((A + 1) + (A - 1) * this.cosineOmega + 2 * Math.sqrt(A) * alpha)) / this.a0;
        this.b1 = (-2 * A * ((A - 1) + (A + 1) * this.cosineOmega)) / this.a0;
        this.b2 = (A * ((A + 1) + (A - 1) * this.cosineOmega - 2 * Math.sqrt(A) * alpha)) / this.a0;
      }
    }

    get coef()
    {
      this.a0 = 1.0;
      return {numerator:[this.b0, this.b1, this.b2], denominator:[this.a0, this.a1, this.a2]};
    }
}

class biquad
{
    constructor(options)
    {
        this._b0 = options.b0;
        this._b1 = options.b1;
        this._b2 = options.b2;

        this._a1 = options.a1;
        this._a2 = options.a2;

        this._prev_input_1 = 0.0;
        this._prev_input_2 = 0.0;
        this._prev_output_1 = 0.0;
        this._prev_output_2 = 0.0;
    }

    process_sample(input_data)
    {
    	let output_data = (this._b0 * input_data          ) + 
                        (this._b1 * this._prev_input_1  ) + 
                        (this._b2 * this._prev_input_2  ) - 
                        (this._a1 * this._prev_output_1 ) - 
                        (this._a2 * this._prev_output_2 );

	    this._prev_input_2  = this._prev_input_1;
	    this._prev_input_1  = input_data;
	    this._prev_output_2 = this._prev_output_1;
	    this._prev_output_1 = output_data;
      return output_data;
    }
}

class WindSourceProcessor extends AudioWorkletProcessor 
{
  constructor(options)
  {
    super();

    this.sampleRate = sampleRate;

    this._biquad = new coefficient({filterType: "band-pass",
                                    frequencyCut: 333,
                                    sampleRate: this.sampleRate,
                                    res: 1.414,
                                    slope: 0,
                                    magnitude: 0});

    this._bp333 = [];
    this._bp333.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));
    this._bp333.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));                                   

    this._biquad.filterType = 'low-pass1';
    this._biquad.frequencyCut = 16;        
    this._lp16 = [];
    this._lp16.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));
    this._lp16.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));   

    this._biquad.filterType = 'low-pass1';
    this._biquad.frequencyCut = 12000;        
    this._lp12000 = [];
    this._lp12000.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));
    this._lp12000.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));   
  }

  process (inputs, outputs, parameters) 
  {
    const output = outputs[0];
    output.forEach((channel, index) => 
    {
      for (let i = 0; i < channel.length; i++) 
      {
        let noise0 = gaussian();
        let noise1 = gaussian();
        let temp1 = this._lp16[index].process_sample(noise0);
        channel[i] = this._lp12000[index].process_sample(this._bp333[index].process_sample(temp1 * noise1) * 0.025);
      }
    })
    return true;
  }
}
  
registerProcessor('wind-source-processor', WindSourceProcessor);

class DripSourceProcessor extends AudioWorkletProcessor
{
  constructor(options)
  {
    super();

    this.sampleRate = sampleRate;

    this._biquad = new coefficient({filterType: "low-pass1",
                                    frequencyCut: 16,
                                    sampleRate: this.sampleRate,
                                    res: 0.707,
                                    slope: 0,
                                    magnitude: 0});                            


    this._lp16 = [];
    this._lp16.push(new biquad({ b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));
    this._lp16.push(new biquad({ b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));
        
    this._biquad.frequencyCut = 0.01;
    this._lp001 = [];
    this._lp001.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));
    this._lp001.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));
    
    this._biquad.frequencyCut = 1900;
    this._lp1900 = [];
    this._lp1900.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));
    this._lp1900.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));
        
    this._biquad.filterType = 'high-pass1';
    this._biquad.frequencyCut = 10000;   
    this._hp10000 = [];
    this._hp10000.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));
    this._hp10000.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));


    this._biquad.filterType = 'low-pass2';
    this._biquad.frequencyCut = 3333;
    this._biquad.res = 30;
    this._vcf3333 = [];
    this._vcf3333.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));
    this._vcf3333.push(new biquad({b0: this._biquad.b0, 
                      b1: this._biquad.b1, 
                      b2: this._biquad.b2, 
                      a1: this._biquad.a1, 
                      a2: this._biquad.a2}));

  }

  process (inputs, outputs, parameters) 
  {
    const output = outputs[0];
    output.forEach((channel, index) => 
    {
      for (let i = 0; i < channel.length; i++) 
      {
        let noise1 = gaussian();
        let noise2 = gaussian();
        let noise3 = gaussian();
        let temp1 = this._lp16[index].process_sample(noise1);
        let temp2 = this._hp10000[index].process_sample(this._lp001[index].process_sample(noise2) / this._lp1900[index].process_sample(noise3)) * 0.03;
        this._biquad.frequencyCut = clip((temp1 * 20000) + 1000, 800, 5000);
        this._vcf3333[index]._b0 = this._biquad.b0;
        this._vcf3333[index]._b1 = this._biquad.b1;
        this._vcf3333[index]._b2 = this._biquad.b2;
        this._vcf3333[index]._a1 = this._biquad.a1;
        this._vcf3333[index]._a2 = this._biquad.a2;
        let temp3 = this._vcf3333[index].process_sample(clip(temp2 * temp2, 0, 0.9));
        channel[i] = temp3;
      }
    })
    return true;
  }
}

registerProcessor('drip-source-processor', DripSourceProcessor);
  
class DropSourceProcessor extends AudioWorkletProcessor 
{

	static get parameterDescriptors() 
	{
		return [
			["lowpassCutOff", 7000, 7000, 11000, "k-rate"],	
		].map(x => new Object({
			name:           x[0],
			defaultValue:   x[1],
			minValue:       x[2],
			maxValue:       x[3],
			automationRate: x[4]
		}));
	}

  constructor(options)
  {
    super();
    this.sampleRate = sampleRate;
    this._invSampleRate = 1 / this.sampleRate; 

    this._cycle = [Math.random() * 2, Math.random() * 2];
    this._degree = [0, 0];

    let model = Math.floor(Math.random() * 4 + 1);
    if(model == 1)
    {
      this._ratio = 344;
      this._offset = 3000;
      this._energy = 2.9;
    }
    else if(model == 2)
    {
      this._ratio = 3000;
      this._offset = 2000;
      this._energy = 2.9;
    }
    else
    {
      this._ratio = 5000;
      this._offset = 2000;
      this._energy = 2.9;
    }

    this._ratio = 3000;
    this._offset = 2000;
    this._energy = 2.9;

    this._biquad = new coefficient({filterType: "band-pass",
                                    frequencyCut: 2.064,
                                    sampleRate: this.sampleRate,
                                    res: 0.002,
                                    slope: 0,
                                    magnitude: 0});  

    this._bp0 = [];
    this._bp0.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));
    this._bp0.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));
    
    this._bp1 = []; 
    this._bp1.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));
    this._bp1.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));

    this._biquad.filterType = "high-pass1";
    this._biquad.frequencyCut = 500;
    this._hp0 = []; 
    this._hp0.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));
    this._hp0.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));

    this._biquad.filterType = "high-pass1";
    this._biquad.frequencyCut = 500;
    this._hp1 = []; 
    this._hp1.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));
    this._hp1.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));

    this._biquad.filterType = "low-pass1";
    this._biquad.frequencyCut = 7000;    
    this._lp = []; 
    this._lp.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));
    this._lp.push(new biquad({b0: this._biquad.b0, 
                    b1: this._biquad.b1, 
                    b2: this._biquad.b2, 
                    a1: this._biquad.a1, 
                    a2: this._biquad.a2}));
    }

    process (inputs, outputs, parameters) 
    {
      const output = outputs[0];       
      output.forEach((channel, index) => 
      {
        for (let i = 0; i < channel.length; i++) 
        {
          let temp = this._bp1[index].process_sample(this._bp0[index].process_sample(Math.random() * 2 - 1) * this._energy);
          this._cycle[index] = this._cycle[index] + (temp * this._ratio + this._offset) * this._invSampleRate;
          this._cycle[index] = this._cycle[index] - Math.floor(this._cycle[index]);
          this._degree[index] = this._cycle[index] * Math.PI * 2;
          let oscillator = Math.sin(this._degree[index]) * Math.pow(temp, 8);
          channel[i] = this._lp[index].process_sample(this._hp1[index].process_sample(this._hp0[index].process_sample(oscillator * 0.1)) * 0.01);
        }
      })
      return true;
    }
}
  
registerProcessor('drop-source-processor', DropSourceProcessor);

class GaussianNoiseProcessor extends AudioWorkletProcessor 
{
  constructor(options)
  {
    super();
    this.sampleRate = sampleRate;
    this._biquad = new coefficient({filterType: "low-pass1",
                                    frequencyCut: 11000,
                                    sampleRate: this.sampleRate,
                                    res: 0.707,
                                    slope: 0,
                                    magnitude: 0});  

    this._lp = [];

    this._lp.push(new biquad({b0: this._biquad.b0, 
                  b1: this._biquad.b1, 
                  b2: this._biquad.b2, 
                  a1: this._biquad.a1, 
                  a2: this._biquad.a2}));

    this._lp.push(new biquad({b0: this._biquad.b0, 
                  b1: this._biquad.b1, 
                  b2: this._biquad.b2, 
                  a1: this._biquad.a1, 
                  a2: this._biquad.a2}));
  }

  process (inputs, outputs, parameters) 
  {
    const output = outputs[0];        
    output.forEach((channel, index) => 
    {
      for (let i = 0; i < channel.length; i++) 
      {
        channel[i] = this._lp[index].process_sample(gaussian() * 0.4);
      }
    })
    return true;
  }
}

registerProcessor('gaussian-noise-processor', GaussianNoiseProcessor);

class DattorroReverbProcessor extends AudioWorkletProcessor 
{
	
	static get parameterDescriptors() 
	{
		return [
			["preDelay",        0.00,   0,  sampleRate - 1, "k-rate"],	
			["bandwidth",       0.50,   0,  1,              "k-rate"],	
			["inputDiffusion1", 0.80,   0,  1,              "k-rate"],	
			["inputDiffusion2", 0.80,   0,  1,              "k-rate"],	
			["decay",           0.85,   0,  1,              "k-rate"],	
			["decayDiffusion1", 0.10,   0,  0.999999,       "k-rate"],	
			["decayDiffusion2", 0.10,   0,  0.999999,       "k-rate"],	
			["damping",         0.32,   0,  1,              "k-rate"],	
			["excursionRate",   8.00,   0,  10,             "k-rate"],	
			["excursionDepth",  1.00,   0,  10,             "k-rate"],	
			["wet",             0.42,   0,  1,              "k-rate"],	
			["dry",             0.58,   0,  1,              "k-rate"],
      ["mono",            0,      0,  1,              "k-rate"],
		].map(x => new Object({
			name:           x[0],
			defaultValue:   x[1],
			minValue:       x[2],
			maxValue:       x[3],
			automationRate: x[4]
		}));
	}

	constructor(options) 
	{
		super(options); 
    
		this._Delays    = [];
		this._pDLength  = sampleRate + (128 - sampleRate % 128); // Pre-delay is always one-second long, rounded to the nearest 128-chunk
		this._preDelay  = new Float32Array(this._pDLength);
		this._pDWrite   = 0;
		this._lp1       = 0.0;
		this._lp2       = 0.0;
		this._lp3       = 0.0;
		this._excPhase	= 0.0;

		[
			0.004771345, 0.003595309, 0.012734787, 0.009307483, 
			0.022579886, 0.149625349, 0.060481839, 0.124995855, 
			0.030509727, 0.141695508, 0.089244313, 0.106280031
		].forEach(x => this.makeDelay(x));

		this._taps = Int16Array.from([
			0.008937872, 0.099929438, 0.064278754, 0.067067639, 0.066866033, 0.006283391, 0.035818689, 
			0.011861161, 0.121870905, 0.041262054, 0.08981553 , 0.070931756, 0.011256342, 0.004065724
		], x => Math.round(x * sampleRate));
	}

	makeDelay(length) 
	{ 
		// len, array, write, read, mask
		let len = Math.round(length * sampleRate);
		let nextPow2 = 2**Math.ceil(Math.log2((len)));
		this._Delays.push([
			new Float32Array(nextPow2),
			len - 1,
			0|0,
			nextPow2 - 1
		]);
	}

	writeDelay(index, data) 
	{
		return this._Delays[index][0][this._Delays[index][1]] = data;
	}

	readDelay(index) 
	{
		return this._Delays[index][0][this._Delays[index][2]];
	}

	readDelayAt(index, i) 
	{
		let d = this._Delays[index];
		return d[0][(d[2] + i)&d[3]];
	}

	// cubic interpolation
	// O. Niemitalo: https://www.musicdsp.org/en/latest/Other/49-cubic-interpollation.html
	readDelayCAt(index, i) 
	{ 
		let d    = this._Delays[index],
			  frac = i-~~i,
			  int  = ~~i + d[2] - 1,
			  mask = d[3];

		let x0 = d[0][int++ & mask],
			  x1 = d[0][int++ & mask],
			  x2 = d[0][int++ & mask],
			  x3 = d[0][int   & mask];

		let a  = (3*(x1-x2) - x0 + x3) / 2,
			  b  = 2*x2 + x0 - (5*x1+x3) / 2,
			  c  = (x2-x0) / 2;

		return (((a * frac) + b) * frac + c) * frac + x1;
	}

	// First input will be downmixed to mono if number of channels is not 2
	// Outputs Stereo.
	process(inputs, outputs, parameters) 
	{
		const 	pd   = ~~parameters.preDelay[0]                           ,
				    bw   = parameters.bandwidth[0]                            ,
				    fi   = parameters.inputDiffusion1[0]                      , 
				    si   = parameters.inputDiffusion2[0]                      ,
				    dc   = parameters.decay[0]                                ,
				    ft   = parameters.decayDiffusion1[0]                      ,
				    st   = parameters.decayDiffusion2[0]                      ,
				    dp   = 1 - parameters.damping[0]                          ,
				    ex   = parameters.excursionRate[0]   / sampleRate         ,
				    ed 	 = parameters.excursionDepth[0]  * sampleRate / 1000  ,
				    we   = parameters.wet[0]             * 0.6                , // lo & ro both mult. by 0.6 anyways
				    dr   = parameters.dry[0]                                  ,
            mn   = parameters.mono[0]                                 ;

		// write to predelay and dry output
		if (inputs[0].length == 2) 
		{
			for (let i = 127; i >= 0; i--) 
			{
				this._preDelay[this._pDWrite+i] = (inputs[0][0][i] + inputs[0][1][i]) * 0.5;

        if(mn == 1)
        {
          outputs[0][0][i] = (inputs[0][0][i]*dr + inputs[0][1][i]*dr) * 0.5;
				  outputs[0][1][i] = (inputs[0][0][i]*dr + inputs[0][1][i]*dr) * 0.5;
        }
        else
        {
          outputs[0][0][i] = inputs[0][0][i]*dr;
				  outputs[0][1][i] = inputs[0][1][i]*dr;
        }
			}
		} 
		else if (inputs[0].length > 0) 
		{
			this._preDelay.set(
				inputs[0][0],
				this._pDWrite
			);
			for (let i = 127; i >= 0; i--) 
				outputs[0][0][i] = outputs[0][1][i] = inputs[0][0][i]*dr;
		} 
		else 
		{
			this._preDelay.set(
        new Float32Array(128),
				this._pDWrite
			);
		}

		let i = 0|0;
		while (i < 128) 
		{
			let lo = 0.0,
				  ro = 0.0;

			this._lp1 += bw * (this._preDelay[(this._pDLength + this._pDWrite - pd + i)%this._pDLength] - this._lp1);

			// pre-tank
			let pre = this.writeDelay(0,             this._lp1          - fi * this.readDelay(0) );
				  pre = this.writeDelay(1, fi * (pre - this.readDelay(1)) +      this.readDelay(0) );
				  pre = this.writeDelay(2, fi *  pre + this.readDelay(1)  - si * this.readDelay(2) );
				  pre = this.writeDelay(3, si * (pre - this.readDelay(3)) +      this.readDelay(2) );

			let split = si * pre + this.readDelay(3);

			// excursions
			let exc   = ed * (1 + Math.cos(this._excPhase * 2 * Math.PI)); 
			let exc2  = ed * (1 + Math.sin(this._excPhase * 2 * Math.PI)); 
			
			// left loop
			let temp =  this.writeDelay( 4, split + dc * this.readDelay(11)    + ft * this.readDelayCAt(4, exc) ); // tank diffuse 1
						      this.writeDelay( 5,         this.readDelayCAt(4, exc)  - ft * temp                      ); // long delay 1
						      this._lp2      += dp * (this.readDelay(5) - this._lp2)                                   ; // damp 1
				  temp =  this.writeDelay( 6,         dc * this._lp2             - st * this.readDelay(6)         ); // tank diffuse 2
						      this.writeDelay( 7,         this.readDelay(6)          + st * temp                      ); // long delay 2
			// right loop 
				  temp =  this.writeDelay( 8, split + dc * this.readDelay(7)     + ft * this.readDelayCAt(8, exc2)); // tank diffuse 3
						      this.writeDelay( 9,         this.readDelayCAt(8, exc2) - ft * temp                      ); // long delay 3
						      this._lp3      += dp * (this.readDelay(9) - this._lp3)                                   ; // damp 2
				  temp =	this.writeDelay(10,         dc * this._lp3             - st * this.readDelay(10)        ); // tank diffuse 4
						      this.writeDelay(11,         this.readDelay(10)         + st * temp                      ); // long delay 4

			lo =  this.readDelayAt( 9, this._taps[0])
				  + this.readDelayAt( 9, this._taps[1])
				  - this.readDelayAt(10, this._taps[2])
				  + this.readDelayAt(11, this._taps[3])
				  - this.readDelayAt( 5, this._taps[4])
				  - this.readDelayAt( 6, this._taps[5])
				  - this.readDelayAt( 7, this._taps[6]);

			ro =  this.readDelayAt( 5, this._taps[7])
				  + this.readDelayAt( 5, this._taps[8])
				  - this.readDelayAt( 6, this._taps[9])
				  + this.readDelayAt( 7, this._taps[10])
				  - this.readDelayAt( 9, this._taps[11])
				  - this.readDelayAt(10, this._taps[12])
				  - this.readDelayAt(11, this._taps[13]);

			outputs[0][0][i] += lo * we;
			outputs[0][1][i] += ro * we;		

			this._excPhase  += ex;	

			i++;

			for (let j = 0, d = this._Delays[0]; j < this._Delays.length; d = this._Delays[++j]) 
			{
				d[1] = (d[1] + 1) & d[3];
				d[2] = (d[2] + 1) & d[3]; 
			}
		}

		// Update preDelay index
		this._pDWrite = (this._pDWrite + 128) % this._pDLength;

		return true;
	}
}

registerProcessor('dattorro-reverb-processor', DattorroReverbProcessor);

class PanKahProcessor extends AudioWorkletProcessor 
{
  static get parameterDescriptors() 
	{
		return [
			["MOD_A",       0.75, 0.00, 1.00, "k-rate"],	
			["MOD_B",       0.87, 0.00, 1.00, "k-rate"],	
			["PITCH_A",     0.55, 0.00, 1.00, "k-rate"],	
			["PITCH_B",     0.95, 0.00, 1.00, "k-rate"],	
			["FREQ_F",      0.89, 0.00, 1.00, "k-rate"],	
			["RES_F",       0.75, 0.00, 1.00, "k-rate"],	
			["CHAOS_A",     0.96, 0.00, 1.00, "k-rate"],	
			["CHAOS_B",     0.31, 0.00, 1.00, "k-rate"],	
			["CHAOS_F",     0.77, 0.00, 1.00, "k-rate"],	
			["MODE_F",      1,    1,    3,    "k-rate"],	
			["MODE_OUT",    5,    1,    5,    "k-rate"],	
			["VOLUME",      0.10, 0.00, 1.00, "k-rate"],	
      
		].map(x => new Object({
			name:           x[0],
			defaultValue:   x[1],
			minValue:       x[2],
			maxValue:       x[3],
			automationRate: x[4]
		}));
	}

  constructor(options)
  {
    super();
    this.sampleRate = sampleRate;
    this.osc_a = new functionGenerator(options);
    this.osc_b = new functionGenerator(options);

    this._biquad = new coefficient({filterType: "high-pass2",
                                    frequencyCut: 150,
                                    sampleRate: this.sampleRate,
                                    res: 0.707,
                                    slope: 0,
                                    magnitude: 0});

    this._hp = new biquad({ b0: this._biquad.b0, 
                            b1: this._biquad.b1, 
                            b2: this._biquad.b2, 
                            a1: this._biquad.a1, 
                            a2: this._biquad.a2});      

    this.svf = new svf(options, 1200, 1.414);

    this.osc_a._cycle = Math.random() / 2;
    this.osc_b._cycle = Math.random() / 2;
    this._invSampleRate = 1 / this.sampleRate;

    this._temp_a = null;
    this._temp_b = null;

    this._chaos = new shift_register();
    this._xor = 1; 
    this._clock = 1;

    this._summing_a = 0;
    this._summing_b = 0;
    this._summing_f = 0;

    this._out = null;
  }

  process (inputs, outputs, parameters) 
  {
    const output = outputs[0];        
    output.forEach((channel) => 
    {
      let tolerance = Math.random() / 10000;
      let key1 = parameters.MODE_OUT[0];
      let key2 = parameters.MODE_F[0];
      for (let i = 0; i < channel.length; i++) 
      {
        let chaos_cv  = (this._chaos.buffer[5] * 0.5) + (this._chaos.buffer[6] * 0.25) + (this._chaos.buffer[7] * 0.125);
        let a_cv      = clip(this._summing_a, 0, 9);
        let b_cv      = clip(this._summing_b, 0, 9);
        let filter_cv = clip(this._summing_f, 0, 9);
        this.osc_a._pitch = (vca(a_cv)+tolerance) * 5000.0;
        this.osc_b._pitch = (vca(b_cv)+tolerance) * 5000.0;
        this.svf.cutoffFrequency = (filter_cv+tolerance) * 4800.0;
        this.svf.resonance = Math.pow(4, parameters.RES_F[0]+tolerance);
        this._temp_a = this.osc_a.Generate();
        this._temp_b = this.osc_b.Generate();
        
        this._xor = this._chaos.xor(this._temp_a.square, this._chaos.buffer[7]);
        this._chaos.sample(this._temp_b.square, this._xor);
        this._summing_a = (1 - parameters.PITCH_A[0]) * 3.0 - (parameters.CHAOS_A[0] * chaos_cv * 1.8) - (this._temp_b.triangle * parameters.MOD_A[0] * 1.8);
        this._summing_b = (1 - parameters.PITCH_B[0]) * 3.0 - (parameters.CHAOS_B[0] * chaos_cv * 1.8) - (this._temp_a.triangle * parameters.MOD_B[0] * 1.8);

        this._summing_f = (parameters.FREQ_F[0]) * 1.8 - (parameters.CHAOS_F[0] * chaos_cv * 1.8);
        let pwm = comparator(this._temp_a.triangle, this._temp_b.triangle);

        switch (key1) 
        {
          case 1:
            this.out = this._temp_a.triangle;
            break;
          case 2:
            this.out = this._temp_b.triangle;
            break;
          
          case 3:
            this.out = this._temp_a.square;
            break;
          
          case 4:
            this.out = this._temp_b.square;
            break;
          
          case 5:
            this.out = pwm;
            break;
        
          default:
            this.out = pwm;
            break;
        }

        switch (key2) {
          case 1:
            channel[i] = this._hp.process_sample(this.svf.process(this.out).LP * parameters.VOLUME[0]);
            break;
          
          case 2:
            channel[i] = this._hp.process_sample(this.svf.process(this.out).BP * parameters.VOLUME[0]);
            break;

          case 3:
            channel[i] = this._hp.process_sample(this.svf.process(this.out).HP * parameters.VOLUME[0]);
            break;
        
          default:
            channel[i] = this._hp.process_sample(this.svf.process(this.out).LP * parameters.VOLUME[0]);
            break;
        }
      }
    })
    return true;
  }
}

registerProcessor('pan-kah-processor', PanKahProcessor);

class LevelMeterProcessor extends AudioWorkletProcessor
{

  constructor(options)
  {
    super();
    this.sampleRate = sampleRate;               
  }

  process(inputList, outputList, parameters) 
  {
    let energy = 0;
    let peak = 0;
    let temp = 0;
    const sourceLimit = Math.min(inputList.length, outputList.length);
    for (let inputNum = 0; inputNum < sourceLimit; inputNum++) 
    {
      let input = inputList[inputNum];
      let output = outputList[0];
      let channelCount = Math.min(input.length, output.length);
      for (let channelNum = 0; channelNum < channelCount; channelNum++) 
      {
        let sampleCount = input[channelNum].length;

        for (let i = 0; i < sampleCount; i++) 
        {
          let sample = input[channelNum][i];
          output[channelNum][i] = sample;
          temp = Math.abs(sample);
          if(peak < temp) { peak = temp; }
          energy = energy + temp;
        }
      }
    }
    this.port.postMessage({"energy": (energy * energy) / 256, "peak": peak});
    return true;
  }
}

registerProcessor('level-meter-processor', LevelMeterProcessor);

class KarlplusStrongProcessor extends AudioWorkletProcessor 
{

  constructor(options)
  {
    super(options);
    this.sampleRate = sampleRate;
    this.lpf = new CircularBuffer(Math.ceil(this.sampleRate / mtof(8)));
    this.apf = new CircularBuffer(Math.ceil(this.sampleRate / mtof(8)));
    this.index = 0;
    this.playing = true;

    this.port.onmessage = (e) => 
    {
        if (e.data === 'stop')
            this.playing = false;
    }
  }

  static get parameterDescriptors() 
  {
		return [
      ["time",    0,        0,        20000,        "a-rate"],	
			["tuning",  0,        -10,      10,           "a-rate"],	
			["decay",   0.95,     0.95,     0.99995,      "a-rate"],	
		].map(x => new Object({
      name: x[0],
			defaultValue: x[1],
			minValue: x[2],
			maxValue: x[3],
			automationRate: x[4]
		}));
  } 

  process(inputs, outputs, parameters)
  {
    const input = inputs[0][0];
    const output = outputs[0][0];
        
    for (let i = 0; i < output.length; i++) 
    {
      output[i] = input ? input[i] : 0;
      let period = parameters['time'].length > 1 ? parameters['time'][i] : parameters['time'][0] ;
      output[i] = output[i] + parameters['decay'][0] * (0.5 * this.lpf.read(i + this.index - period) + 0.5 * this.lpf.read(i + this.index - period - 1));
      this.lpf.write(this.index + i, output[i]);

      let temp = parameters['tuning'].length > 1 ? parameters['tuning'][i] : parameters['tuning'][0] ;
      output[i] = temp * output[i] + this.lpf.read(this.index + i - 1) - temp * this.apf.read(this.index + i - 1);

      this.apf.write(this.index + i, output[i]);
    }
    this.index += output.length;
    return this.playing;
  }
}

registerProcessor('karlplus-strong-processor', KarlplusStrongProcessor)

class GrainularProcessor extends AudioWorkletProcessor
{

  constructor(options)
  {
    super(options);
    this.sampleRate = sampleRate;    
    this.fftSize = 1024;
    this.offset = this.fftSize - 128;

    this.inputBuffer = new circularBuffer();
    this.inputBuffer.createBuffer(this.fftSize);
    this.phaseBuffer = new circularBuffer();
    this.phaseBuffer.createBuffer(this.fftSize);

    this.xn = new Array(128).fill(0);

    this.Xn_imag = new Array(this.fftSize).fill(0);
    this.Xn_real = new Array(this.fftSize).fill(0);

    this.x_synthesis_real = new Array(this.fftSize).fill(0);
    this.x_synthesis_imag = new Array(this.fftSize).fill(0);

    this.randomTemp = new Array((this.fftSize / 2) - 1).fill(0);
  }

  createWindow()
  {
    this.window = new Array(this.fftSize).fill(0);
    for (let i = 0; i < this.window.length; i++) 
    {
      let omega = (i - ((this.fftSize - 1) / 2)) / ((this.fftSize - 1) / 2) * Math.PI;
      this.window[i] =  (1 + Math.cos(omega))/2;
    }
  }

  static get parameterDescriptors() 
  {
		return [
      ["bandshift",     0,    0,    255,    "k-rate"],
      ["randomness",    0,    -1,    1,    "k-rate"],
      ["phaseshift",    0.0,   -1,    1,    "k-rate"],
      ["phasenull",     0,    0,    1,    "k-rate"],
      ["panning",       0.5,    0,    1,    "k-rate"],
		].map(x => new Object({
      name: x[0],
			defaultValue: x[1],
			minValue: x[2],
			maxValue: x[3],
			automationRate: x[4]
		}));
  } 

  process(inputs, outputs, parameters)
  {
    if (inputs[0].length > 0)
    {
      let bandshift = parameters['bandshift'][0];
      let randomness = parameters['randomness'][0];
      let phaseshift =  parameters['phaseshift'][0];
      let phasenull =  parameters['phasenull'][0];
      let panning =  parameters['panning'][0];

      for (let i = 0; i < 128; i++)
      {
        let blend = inputs[0][0][i] * panning + inputs[0][1][i] * (1 - panning);
        this.xn[i] = blend;
        this.inputBuffer.writeBuffer(blend);
      }

      // prepare vector for fft
      this.Xn_imag.fill(0);
      for (let i = 0; i < this.fftSize; i++)
      {
        this.Xn_real[i] = this.inputBuffer.readBuffer(this.fftSize - i);
      }

      // perform fft
      fft(this.Xn_real, this.Xn_imag);
      this.Xn_real = complex_roll(this.Xn_real, bandshift);
      // this.Xn_imag = complex_roll(this.Xn_imag, bandshift);

      // prepare mirrored random phase vector
      this.randomTemp.forEach((value, index, array) =>
      {
        array[index] = (Math.random() * 2 * Math.PI - Math.PI);
        // array[index] = gaussian() * Math.PI;
      });
      let flip = this.randomTemp.reverse().map(function(x) { return x * -1; });
      let randomPhase = [0].concat(this.randomTemp).concat([0]).concat(flip);

      // convert fft data into amplitude and angle
      for (let i = 0; i < this.fftSize; i++)
      {
        let amplitude = Math.sqrt(this.Xn_real[i] * this.Xn_real[i] + this.Xn_imag[i] * this.Xn_imag[i]);
        let angle = Math.atan2(this.Xn_imag[i], this.Xn_real[i]) * 360 / (2 * Math.PI);
        angle = (angle + 180 * phaseshift) % 360;
        let angle_hat = ((angle / 360) * 2 * Math.PI + randomPhase[i] * randomness) * (1 - phasenull);

        let [real, imag] = polar2complex(amplitude, angle_hat);
        this.x_synthesis_real[i] = real;
        this.x_synthesis_imag[i] = imag;
      }

      // perform ifft
      ifft(this.x_synthesis_real, this.x_synthesis_imag);

      for (let i = 0; i < 128; i++) 
      {
        outputs[0][0][i] = this.x_synthesis_real[i + this.offset];
      }
      this.port.postMessage({"buffer1": inputs[0][0], "buffer2": inputs[0][1]});
    }
    return true;      
  }
}

registerProcessor('grainular-processor', GrainularProcessor);

class FeedbackDelayNetworkProcessor extends AudioWorkletProcessor
{
  constructor(options)
  {
    super(options);
    this.sampleRate = sampleRate;

    this.wave1 = new functionGenerator();
    this.wave2 = new functionGenerator();
    this.wave3 = new functionGenerator();
    this.wave4 = new functionGenerator();

    this.wave1._cycle = 0;
    this.wave2._cycle = 0.25;
    this.wave3._cycle = 0.5;
    this.wave4._cycle = 0.75;

    this.wave1._pitch = 1;
    this.wave2._pitch = 1;
    this.wave3._pitch = 1;
    this.wave4._pitch = 1;

    this.cb1 = new circularBuffer();
    this.cb2 = new circularBuffer();
    this.cb3 = new circularBuffer();
    this.cb4 = new circularBuffer();

    this.cb1.createBuffer(4096);
    this.cb2.createBuffer(4096);
    this.cb3.createBuffer(4096);
    this.cb4.createBuffer(4096);

    this._coefficient = new coefficient(
      {
        filterType: "low-pass1",
        frequencyCut: 150,
        sampleRate: this.sampleRate,
        res: 0,
        slope: 0,
        magnitude: 0
      });

    this._biquad1 = new biquad(
      {
        b0: this._coefficient.b0, 
        b1: this._coefficient.b1, 
        b2: this._coefficient.b2, 
        a1: this._coefficient.a1, 
        a2: this._coefficient.a2
      });

    this._biquad2 = new biquad(
      {
        b0: this._coefficient.b0, 
        b1: this._coefficient.b1, 
        b2: this._coefficient.b2, 
        a1: this._coefficient.a1, 
        a2: this._coefficient.a2
      });

    this._biquad3 = new biquad(
      {
        b0: this._coefficient.b0, 
        b1: this._coefficient.b1, 
        b2: this._coefficient.b2, 
        a1: this._coefficient.a1, 
        a2: this._coefficient.a2
      });

    this._biquad4 = new biquad(
      {
        b0: this._coefficient.b0, 
        b1: this._coefficient.b1, 
        b2: this._coefficient.b2, 
        a1: this._coefficient.a1, 
        a2: this._coefficient.a2
      });

  }

  static get parameterDescriptors() 
  {
		return [
      ["bypass",        0,      0,    1,    "k-rate"],
      ["brightness",    1200,   150,  5000, "k-rate"],
      ["damping",       0.5,    0,    1,    "k-rate"],
      ["decay",         0.9,    0,    1,    "k-rate"],
      ["size",          0.9,    0,    1,    "k-rate"],
      ["frequency",     2,      0,    20,   "k-rate"],
      ["depth",         63,     0,    127,  "k-rate"],
      ["mix",           0.5,    0,    1,    "k-rate"],
      ["freeze",        0,      0,    1,    "k-rate"],
		].map(x => new Object({
      name: x[0],
			defaultValue: x[1],
			minValue: x[2],
			maxValue: x[3],
			automationRate: x[4]
		}));
  } 

  process(inputs, outputs, parameters)
  {
    if (inputs[0].length == 2)
    {
      let bypass      = parameters['bypass'][0];
      let brightness  = parameters['brightness'][0];
      let damping     = parameters['damping'][0];
      let decay       = parameters['decay'][0];
      let size        = parameters['size'][0];
      let frequency   = parameters['frequency'][0];
      let depth       = parameters['depth'][0];
      let mix         = parameters['mix'][0];
      let freeze      = parameters['freeze'][0];
      if(freeze == 1)
      {

      }
      if(bypass == 1) { mix = 0; }

      this.wave1._pitch = frequency;
      this.wave2._pitch = frequency;
      this.wave3._pitch = frequency;
      this.wave4._pitch = frequency;

      this._coefficient.frequencyCut = brightness;
      this._biquad1._a1 = this._coefficient.a1;
      this._biquad1._a2 = this._coefficient.a2;
      this._biquad1._b0 = this._coefficient.b0;
      this._biquad1._b1 = this._coefficient.b1;
      this._biquad1._b2 = this._coefficient.b2;

      this._biquad2._a1 = this._coefficient.a1;
      this._biquad2._a2 = this._coefficient.a2;
      this._biquad2._b0 = this._coefficient.b0;
      this._biquad2._b1 = this._coefficient.b1;
      this._biquad2._b2 = this._coefficient.b2;

      this._biquad3._a1 = this._coefficient.a1;
      this._biquad3._a2 = this._coefficient.a2;
      this._biquad3._b0 = this._coefficient.b0;
      this._biquad3._b1 = this._coefficient.b1;
      this._biquad3._b2 = this._coefficient.b2;

      this._biquad4._a1 = this._coefficient.a1;
      this._biquad4._a2 = this._coefficient.a2;
      this._biquad4._b0 = this._coefficient.b0;
      this._biquad4._b1 = this._coefficient.b1;
      this._biquad4._b2 = this._coefficient.b2;

      for (let i = 0; i < 128; i++)
      {
        let drySignal1 = inputs[0][0][i] * (1 - bypass);
        let drySignal2 = inputs[0][1][i] * (1 - bypass);

        let modulation1 = this.wave1.GenerateSine() * depth;
        let modulation2 = this.wave2.GenerateSine() * depth;
        let modulation3 = this.wave3.GenerateSine() * depth;
        let modulation4 = this.wave4.GenerateSine() * depth;

        let feedback1 = this.cb1.readBuffer((2819.0 + modulation1) * size);
        let feedback2 = this.cb2.readBuffer((3343.0 + modulation2) * size);
        let feedback3 = this.cb3.readBuffer((3581.0 + modulation3) * size);
        let feedback4 = this.cb4.readBuffer((4133.0 + modulation4) * size);

        let lpf1 = this._biquad1.process_sample(feedback1);
        let lpf2 = this._biquad2.process_sample(feedback2);
        let lpf3 = this._biquad3.process_sample(feedback3);
        let lpf4 = this._biquad4.process_sample(feedback4);

        let damp1 = (lpf1 - feedback1) * damping;
        let damp2 = (lpf2 - feedback2) * damping;
        let damp3 = (lpf3 - feedback3) * damping;
        let damp4 = (lpf4 - feedback4) * damping;

        let A = (damp1 + feedback1) * 0.5 * (decay * 0.25 + 0.75) + drySignal1;
        let B = (damp2 + feedback2) * 0.5 * (decay * 0.25 + 0.75) + drySignal2;
        let C = (damp3 + feedback3) * 0.5 * (decay * 0.25 + 0.75);
        let D = (damp4 + feedback4) * 0.5 * (decay * 0.25 + 0.75);

        let output1 = A + B + C + D;
        let output2 = A - B + C - D;
        let output3 = A + B - C - D;
        let output4 = A - B - C + D;

        this.cb1.writeBuffer(output1);
        this.cb2.writeBuffer(output2);
        this.cb3.writeBuffer(output3);
        this.cb4.writeBuffer(output4);

        outputs[0][0][i] = C * mix + inputs[0][0][i] * (1 - mix);
        outputs[0][1][i] = D * mix + inputs[0][1][i] * (1 - mix);
      }
    }
    return true;
  }
}

registerProcessor('feedback-delay-network-processor', FeedbackDelayNetworkProcessor);

class JohnStautnerMillerPucketteFeedbackDelayNetworkProcessor extends AudioWorkletProcessor
{
  constructor(options)
  {
    super(options);
    this.sampleRate = sampleRate;

    this.wave1 = new functionGenerator();
    this.wave2 = new functionGenerator();
    this.wave3 = new functionGenerator();
    this.wave4 = new functionGenerator();

    this.wave1._cycle = 0;
    this.wave2._cycle = 0.25;
    this.wave3._cycle = 0.5;
    this.wave4._cycle = 0.75;

    this.wave1._pitch = 1;
    this.wave2._pitch = 1;
    this.wave3._pitch = 1;
    this.wave4._pitch = 1;

    this.cb1 = new circularBuffer();
    this.cb2 = new circularBuffer();
    this.cb3 = new circularBuffer();
    this.cb4 = new circularBuffer();

    this.cb1.createBuffer(1024);
    this.cb2.createBuffer(2048);
    this.cb3.createBuffer(4096);
    this.cb4.createBuffer(8192);

    this._coefficient = new coefficient(
      {
        filterType: "low-pass2",
        frequencyCut: 150,
        sampleRate: this.sampleRate,
        res: 0.25,
        slope: 0,
        magnitude: 0
      });

    this._biquad1 = new biquad(
      {
        b0: this._coefficient.b0, 
        b1: this._coefficient.b1, 
        b2: this._coefficient.b2, 
        a1: this._coefficient.a1, 
        a2: this._coefficient.a2
      });

    this._biquad2 = new biquad(
      {
        b0: this._coefficient.b0, 
        b1: this._coefficient.b1, 
        b2: this._coefficient.b2, 
        a1: this._coefficient.a1, 
        a2: this._coefficient.a2
      });
  }

  static get parameterDescriptors() 
  {
		return [
      ["bypass",        0,      0,    1,    "k-rate"],
      ["brightness",    2500,   150,  5000, "k-rate"],
      ["spread",        0.25,   0,    1,    "k-rate"],
      ["decay",         0.91,   0,    1,    "k-rate"],
      ["frequency",     0.375,  0,    1,    "k-rate"],
      ["depth",         127,    0,    127,  "k-rate"],
      ["mix",           1,      0,    1,    "k-rate"],
      ["freeze",        0,      0,    1,    "k-rate"],
		].map(x => new Object({
      name: x[0],
			defaultValue: x[1],
			minValue: x[2],
			maxValue: x[3],
			automationRate: x[4]
		}));
  } 

  process(inputs, outputs, parameters)
  {
    let bypass      = parameters['bypass'][0];
    let brightness  = parameters['brightness'][0];
    let spread      = parameters['spread'][0];
    let decay       = parameters['decay'][0];
    let frequency   = parameters['frequency'][0];
    let depth       = parameters['depth'][0];
    let mix         = parameters['mix'][0];
    let freeze      = parameters['freeze'][0];

    let factor = 2.5;
    this.wave1._pitch = frequency / factor;
    this.wave2._pitch = frequency;
    this.wave3._pitch = frequency * factor;
    this.wave4._pitch = frequency * factor * factor;

    let rad_0 = 0.1250 * 2 * Math.PI;
    let rad_1 = spread * 2 * Math.PI;
    this._coefficient.frequencyCut = brightness;
    this._biquad1._a1 = this._coefficient.a1;
    this._biquad1._a2 = this._coefficient.a2;
    this._biquad1._b0 = this._coefficient.b0;
    this._biquad1._b1 = this._coefficient.b1;
    this._biquad1._b2 = this._coefficient.b2;

    this._biquad2._a1 = this._coefficient.a1;
    this._biquad2._a2 = this._coefficient.a2;
    this._biquad2._b0 = this._coefficient.b0;
    this._biquad2._b1 = this._coefficient.b1;
    this._biquad2._b2 = this._coefficient.b2;

    if (inputs[0].length == 2)
    {
      for (let i = 0; i < 128; i++)
      {
        let drySignal1 = inputs[0][0][i] * (1 - bypass);
        let drySignal2 = inputs[0][1][i] * (1 - bypass);

        let modulation1 = this.wave1.Generate().triangle * depth;
        let modulation2 = this.wave2.Generate().triangle * depth;
        let modulation3 = this.wave3.Generate().triangle * depth;
        let modulation4 = this.wave4.Generate().triangle * depth;
        
        let feedback1 = this.cb1.readBufferLinear((887.0 + modulation1));
        let feedback2 = this.cb2.readBufferLinear((1913.0 + modulation2));
        let feedback3 = this.cb3.readBufferLinear((3967.0 + modulation3));
        let feedback4 = this.cb4.readBufferLinear((8053.0 + modulation4));

        let temp1 = feedback1 * Math.cos(rad_0) - feedback2 * Math.sin(rad_0);
        let temp2 = feedback1 * Math.sin(rad_0) + feedback2 * Math.cos(rad_0);
        let temp3 = feedback3 * Math.cos(rad_0) - feedback4 * Math.sin(rad_0);
        let temp4 = feedback3 * Math.sin(rad_0) + feedback4 * Math.cos(rad_0);
        
        let output1 = temp1 * Math.cos(rad_1) - temp3 * Math.sin(rad_1);
        let output2 = temp1 * Math.sin(rad_1) + temp3 * Math.cos(rad_1);
        let output3 = temp2 * Math.cos(rad_1) - temp4 * Math.sin(rad_1);
        let output4 = temp2 * Math.sin(rad_1) + temp4 * Math.cos(rad_1);

        this.cb1.writeBuffer(drySignal1 + output1 * decay);
        this.cb2.writeBuffer(drySignal2 + output3 * decay);
        this.cb3.writeBuffer(this._biquad1.process_sample(output2 * decay));
        this.cb4.writeBuffer(this._biquad2.process_sample(output4 * decay));

        outputs[0][0][i] = (output3 * 2) * mix + inputs[0][0][i] * (1 - mix);
        outputs[0][1][i] = (output4 * 2) * mix + inputs[0][1][i] * (1 - mix);
      }
    }
    return true;
  }
}

registerProcessor('js-mp-feedback-delay-network-processor', JohnStautnerMillerPucketteFeedbackDelayNetworkProcessor);