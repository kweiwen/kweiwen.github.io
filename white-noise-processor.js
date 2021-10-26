class WhiteNoiseProcessor extends AudioWorkletProcessor {
    constructor(options)
    {
        super();
        this.sampleRate = options.processorOptions.sampleRate;
        this.channelCount = options.processorOptions.channelCount;      
    }

    process (inputs, outputs, parameters) 
    {
        const output = outputs[0];        
        output.forEach((channel) => 
        {
            for (let i = 0; i < channel.length; i++) 
            {
                channel[i] = (Math.random() * 2 - 1);
            }
        })
        return true
    }
}
  
registerProcessor('WhiteNoiseProcessor', WhiteNoiseProcessor)