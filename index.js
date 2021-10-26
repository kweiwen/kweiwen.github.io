let DEBUG_MAIN = false;
let clickTimes = 0;

class BehaviorManager
{
  constructor(options)
  {
    
  }
  
  async init()
  {
    this.context = new AudioContext(); 
    await this.context.audioWorklet.addModule('./drop-processor.js');
    
    let channelCount = this.context.destination.channelCount;
    let sampleRate = this.context.sampleRate;
    let argument = { processorOptions: { channelCount, sampleRate }};
    this.dropNode = new AudioWorkletNode(this.context, 'DropProcessor', argument);
    this.dropNode.connect(this.context.destination);
  }
}

const bm = new BehaviorManager();

document.getElementById("toggle").addEventListener("click", function()
{
  let isChecked = document.getElementById("toggle").checked;
  clickTimes = clickTimes + 1;

  if(clickTimes == 1 & isChecked)
  {
    bm.init();
  }
  else if(isChecked)
  {
    bm.context.resume().then(function() 
    {
      console.log(bm.context.state);
    });
  }
  else
  {
    bm.context.suspend().then(function() 
    {
      console.log(bm.context.state);
    });
  }
});

// let DEBUG_MAIN = false;

// async function createReverb(audioContext, path) 
// {
//     let convolver = audioContext.createConvolver();

//     let response     = await fetch(path);
//     let arraybuffer  = await response.arrayBuffer();
//     convolver.buffer = await audioContext.decodeAudioData(arraybuffer);

//     return convolver;
// }

// async function main() 
// {
//   const audioContext = new AudioContext(); 
//   await audioContext.audioWorklet.addModule('./drop-processor.js');
//   let reverb1 = await createReverb(audioContext, "./IMreverbs/Cement Blocks 1.wav");
//   let reverb2 = await createReverb(audioContext, "./IMreverbs/Cement Blocks 2.wav");
//   reverb1.connect(reverb2);
//   reverb2.connect(audioContext.destination);

//   if(DEBUG_MAIN) 
//   {
//     console.log("Initialization")
//     console.log("Source: ", audioContext);
//   }

//   channelCount = audioContext.destination.channelCount;
//   sampleRate = audioContext.sampleRate;
//   let argument = { processorOptions: { channelCount, sampleRate }}

//   const DropNode1 = new AudioWorkletNode(audioContext, 'DropProcessor', argument);
//   const DropNode2 = new AudioWorkletNode(audioContext, 'DropProcessor', argument);
//   const DropNode3 = new AudioWorkletNode(audioContext, 'DropProcessor', argument);
//   const DropNode4 = new AudioWorkletNode(audioContext, 'DropProcessor', argument);
//   const DropNode5 = new AudioWorkletNode(audioContext, 'DropProcessor', argument);
//   const DropNode6 = new AudioWorkletNode(audioContext, 'DropProcessor', argument);


//   document.getElementById("toggle").addEventListener("click", function()
//   {
//     let isChecked = document.getElementById("toggle").checked;
//     connect(isChecked, reverb1, DropNode1);
//     connect(isChecked, reverb1, DropNode2);
//     connect(isChecked, reverb1, DropNode3);
//     connect(isChecked, reverb1, DropNode4);
//     connect(isChecked, reverb1, DropNode5);
//     connect(isChecked, reverb1, DropNode6);
//   });
// }

// async function connect(isChecked, context, source) 
// {
//   if(DEBUG_MAIN) 
//   {
//     console.log("isChecked: ", isChecked);
//     console.log("Context: ", context);
//     console.log("Source: ", source);
//   }

//   try
//   {
//     if(isChecked)
//       source.connect(context);
//     else
//       source.disconnect(context);
//   }
//   catch(error)
//   {
//     console.error(error);
//   }
// }

// main();