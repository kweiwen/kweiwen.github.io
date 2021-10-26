let DEBUG_SRC_PROCESSOR = true;

class coefficient
{
    constructor(options)
    {

    }

    generate()
    {
        return;
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
    	let output_data = (this._b0 * input_data) + 
        (this._b1 * this._prev_input_1) + 
        (this._b2 * this._prev_input_2) - 
        (this._a1 * this._prev_output_1) - 
        (this._a2 * this._prev_output_2);

	    this._prev_input_2 = this._prev_input_1;
	    this._prev_input_1 = input_data;
	    this._prev_output_2 = this._prev_output_1;
	    this._prev_output_1 = output_data;
        return output_data;
    }
}

class DropProcessor extends AudioWorkletProcessor {
    constructor(options)
    {
        super();
        this.sampleRate = options.processorOptions.sampleRate;
        this.channelCount = options.processorOptions.channelCount;
        this._cycle = 0;
        this._degree = 0;

        this._ratio = 4000 + 400 * (Math.random() * 2 - 1);
        this._offset = 3000 + 300 * (Math.random() * 2 - 1);
        this._energy = 3.3 + 0.33 * (Math.random() * 2 - 1);

        if(DEBUG_SRC_PROCESSOR)
        {
            console.log("ratio: ", this._ratio);
            console.log("offset: ", this._offset);
            console.log("energy: ", this._energy);
        }

        // [0.013328791622752871, 0.0, -0.013328791622752871]
        // [1.0, -1.9733423447318437, 0.9733424167544943]
        this._bp = new biquad({ b0: 0.013328791622752871, 
                                b1: 0.0, 
                                b2: -0.013328791622752871, 
                                a1: -1.9733423447318437, 
                                a2: 0.9733424167544943});

        // [0.23393275413302522, -0.46786550826605044, 0.23393275413302522]
        // [1.0, -0.4673641035094105, -0.5316330869773097]
        this._hp = new biquad({ b0: 0.23393275413302522, 
                                b1: -0.46786550826605044, 
                                b2: 0.23393275413302522, 
                                a1: -0.4673641035094105, 
                                a2: -0.5316330869773097});
        
        // [0.6000028501589929, 0.0, 0.0]
        // [1.0, -0.39999714984100715, 0.0]
        this._lp = new biquad({ b0: 0.6000028501589929, 
                                b1: 0.0, 
                                b2: 0.0, 
                                a1: -0.39999714984100715, 
                                a2: 0.0});

    }

    process (inputs, outputs, parameters) 
    {
        const output = outputs[0];        
        output.forEach((channel) => 
        {
            for (let i = 0; i < channel.length; i++) 
            {
                let background = (Math.random() * 2 - 1) * 0.0005;
                let temp = this._bp.process_sample(Math.random() * 2 - 1) * this._energy;
                this._cycle = this._cycle + (temp * this._ratio + this._offset) / this.sampleRate;
                this._degree = this._cycle * Math.PI * 2;
                let oscillator = Math.sin(this._degree) * Math.pow(temp, 8);
                channel[i] = this._lp.process_sample(this._hp.process_sample(oscillator) + background);
            }
        })
        return true
    }
}
  
registerProcessor('DropProcessor', DropProcessor)
  