class NodeManager
{
  constructor(options)
  {
    this.toggleTime = 0;
    this.context = null;

    this.isChecked = null;
    this.isReady = false;
    this.initNodes();
  }

  initNodes()
  {
    // place audio worklet node here!
    this.trackNode1 = null;
    this.trackNode2 = null;
    this.grainNode = null;
    this.reverbNode  = null;
  }

  async initAudioPreference()
  {
    let audioContext = new AudioContext(); 
    try 
    {
      await audioContext.audioWorklet.addModule("digital-processor.js");
      
      fetch("123456.wav")
      .then(response => response.arrayBuffer())
      .then(arrayBuffer => audioContext.decodeAudioData(arrayBuffer))
      .then(audioBuffer => 
      {
        this.createTopology(audioBuffer, audioContext);
      })
    } 
    catch (e) 
    {
      console.log(e);
    }
    this.context = audioContext;
  }

  createTopology(audioBuffer, audioContext) 
  {    
    this.trackNode1 = audioContext.createBufferSource();
    // this.trackNode2 = audioContext.createBufferSource();

    this.trackNode1.buffer = audioBuffer;
    // this.trackNode2.buffer = audioBuffer;

    this.trackNode1.loop = true;
    // this.trackNode1.loopStart = 1.2531;
    // this.trackNode1.loopEnd = 1.673;
    // this.trackNode1.playbackRate.value = 0.999;

    // this.trackNode2.loop = true;
    // this.trackNode2.loopStart = 1.2531;
    // this.trackNode2.loopEnd = 1.673;
    // this.trackNode2.playbackRate.value = 1.001;

    this.grainNode = new AudioWorkletNode(audioContext, 'grainular-processor', { inputChannelCount: [2], outputChannelCount: [1] });
    this.reverbNode = new AudioWorkletNode(audioContext, 'dattorro-reverb-processor', { inputChannelCount: [2], outputChannelCount: [2] });

    let filter = audioContext.createBiquadFilter();
    filter.type = 'lowpass';
    filter.frequency.value = 5000;

    this.trackNode1.connect(this.grainNode).connect(filter).connect(this.reverbNode).connect(audioContext.destination);
    // this.trackNode2.connect(audioContext.destination);

    this.trackNode1.start();
    // this.trackNode2.start();

    this.isReady = true;
  }

  get MIX() {return this._MIX;}
  set MIX(input_data)
  {
    this._MIX = input_data;
    if(this.context)
    {
      let parameter = this.grainNode.parameters.get("mix");
      parameter.value = input_data / 100;
    }
  }

  get DOMAIN() {return this._DOMAIN;}
  set DOMAIN(input_data)
  {
    this._DOMAIN = input_data;
    if(this.context)
    {
      let parameter = this.grainNode.parameters.get("domain");
      parameter.value = input_data / 100;
    }
  }

  get RANDOMNESS() {return this._RANDOMNESS;}
  set RANDOMNESS(input_data)
  {
    this._RANDOMNESS = input_data;
    if(this.context)
    {
      let parameter = this.grainNode.parameters.get("randomness");
      parameter.value = input_data / 100;
    }
  }

  get START() {return this._START;}
  set START(input_data)
  {
    this._START = input_data;
    if(this.context)
    {
      this.trackNode1.loopStart = (input_data / 100) * this.trackNode1.buffer.duration;
      this.trackNode2.loopStart = (input_data / 100) * this.trackNode2.buffer.duration;
    }
  }

  get LENGTH() {return this._PANNING;}
  set LENGTH(input_data)
  {
    this._LENGTH = input_data;
    if(this.context)
    {
      this.trackNode1.loopEnd = (input_data / 100) * this.trackNode1.buffer.duration;
      this.trackNode2.loopEnd = (input_data / 100) * this.trackNode2.buffer.duration;
    }
  }

  get PANNING() {return this._PITCH;}
  set PANNING(input_data)
  {
    this._PANNING = input_data;
    if(this.context)
    {
      let parameter = this.grainNode.parameters.get("panning");
      parameter.value = input_data / 100;
    }
  }

  get PITCH() {return this._PITCH;}
  set PITCH(input_data)
  {
    this._PITCH = input_data;
    if(this.context)
    {
      let frequency1 = Math.pow(2, (input_data / 12.0));
      let frequency2 = Math.pow(2, (input_data * 1.001 / 12.0));
      this.trackNode1.playbackRate.value = frequency1;
      this.trackNode2.playbackRate.value = frequency2;
    }
  }

  get PHASESHIFT() {return this._PHASESHIFT;}
  set PHASESHIFT(input_data)
  {
    this._PHASESHIFT = input_data;
    if(this.context)
    {
      let parameter = this.grainNode.parameters.get("phaseshift");
      parameter.value = input_data / 100;
    }
  }


  get PHASENULL() {return this._PHASENULL;}
  set PHASENULL(input_data)
  {
    this._PHASENULL = input_data;
    if(this.context)
    {
      let parameter = this.grainNode.parameters.get("phasenull");
      parameter.value = input_data / 100;
    }
  }


  get BANDSHIFT() {return this._BANDSHIFT;}
  set BANDSHIFT(input_data)
  {
    this._BANDSHIFT = input_data;
    if(this.context)
    {
      let parameter = this.grainNode.parameters.get("bandshift");
      parameter.value = parseInt(input_data);
    }
  }
}

const STATE = document.getElementById("STATE");
const STATELABEL = document.getElementById("STATE-LABEL");

window.addEventListener('click', onClickEventHandler);

async function onClickEventHandler(event)
{
  await ToggleAction(event);
  // clicked = clicked + 1;
  // if(clicked == 1)
  // {
  //   singleClickTimer = setTimeout(() => { clicked = 0; console.log("single click!"); TactileAction(event); }, 200);
  // } 
  // else if (clicked === 2) 
  // {
  //   clearTimeout(singleClickTimer);
  //   clicked = 0;
  //   console.log("double click!");
  //   await ToggleAction();
  // }
}

async function ToggleAction(event)
{
  if(event.target.id == "STATE-LABEL")
  {
    nodes.toggleTime = nodes.toggleTime + 1;
    nodes.isChecked = nodes.toggleTime % 2;
    if(nodes.toggleTime == 1)
    {
      try
      {
        await nodes.initAudioPreference();
        console.log(nodes.context.state);
        window.setTimeout(function() 
        {
          timeOutHandler = window.setInterval(refreshConsole, 100);
        }, 100);
        STATE.innerHTML = "PAUSE";
      }
      catch(e)
      {
        console.error(e);
      }
    }
    else if(nodes.isChecked)
    {
      nodes.context.resume().then(function() 
      {
        console.log(nodes.context.state);
        window.setTimeout(function() 
        {
          timeOutHandler = window.setInterval(refreshConsole, 100);
        }, 100);

        STATE.innerHTML = "PAUSE";
      });
    }
    else
    {
      nodes.context.suspend().then(function() 
      {
        resetConsole();
        console.log(nodes.context.state);
        clearInterval(timeOutHandler);
        clearInterval(masterClockTimeOut);
      });
    }
  }
}

function resetConsole()
{
  STATE.innerHTML = "PLAY";
  STATELABEL.innerHTML = ">STATE";
}

function refreshConsole()
{
  animationIndex++;
  animationIndex = (animationIndex % 4);
  if(animationIndex == 0)
  {
    STATELABEL.innerHTML = "\\STATE";
  }
  else if(animationIndex == 1)
  {
    STATELABEL.innerHTML = "|STATE";
  }
  else if(animationIndex == 2)
  {
    STATELABEL.innerHTML = "/STATE";
  }
  else if(animationIndex == 3)
  {
    STATELABEL.innerHTML = "-STATE";
  }
}

const nodes = new NodeManager();

let timeOutHandler = null;
let masterClockTimeOut = null;
let tapped = null;
let clicked = null;
let animationIndex = 0;

let SIZE = 128;
let DATA1 = new Array(SIZE).fill(0);
let DATA2 = new Array(SIZE).fill(0);
let SPACE = (2*Math.PI/SIZE);
let FPS = 30;
let OMEGA = (1 / FPS) * 2 * Math.PI;
let DEG = 360;
let ANGLE = 45;
let SPEED = 1;
let AMOUNT = 1;
let RADIUS = 64;
let CG_pixel;
let counter = 0;

let canvas;

function preload() 
{
  CG_pixel = loadFont('CG-pixel-3x5-mono.woff');
}

function setup() 
{
  canvas = createCanvas(windowWidth/4, windowWidth/4, WEBGL);
  // canvas = createCanvas(windowWidth/4, windowWidth/4);
  angleMode(DEGREES);
  frameRate(FPS);

  textFont(CG_pixel);
  textSize(0.05 * windowWidth);
  textAlign(CENTER, CENTER);
}

function draw() 
{
  background(120);
  // background(0);

  stroke(255, 255, 255);
  strokeWeight(6);
  line(windowWidth/8, windowWidth/8, 0, -windowWidth/8, windowWidth/8, 0);
  line(windowWidth/8, -windowWidth/8, 0, -windowWidth/8, -windowWidth/8, 0);
  line(windowWidth/8, windowWidth/8, 0, windowWidth/8, -windowWidth/8, 0);
  line(-windowWidth/8, windowWidth/8, 0, -windowWidth/8, -windowWidth/8, 0);

  if(nodes.isChecked & nodes.isReady)
  {
    gui_disp_fps();
  }


  if(nodes.isChecked & nodes.isReady)
  {
    // gui_disp_circle();
    gui_plot_signal();
  }
}

function windowResized() 
{
  textSize(0.05 * windowWidth);
  resizeCanvas(windowWidth/4, windowWidth/4);
}

function gui_disp_fps()
{
  fill(255, 255, 255, 255);
  text(counter++, 0, 0);
}

function gui_disp_circle()
{
  noFill();
  strokeWeight(1);
  beginShape();
  for(let j = 0; j < 360; j+=10)
  {
    stroke(255, 255, 255);
    let x = RADIUS * cos(j);
    let y = RADIUS * sin(j);
    let z = 0;
    vertex(x, y, z);
  }
  endShape(CLOSE);
}

function gui_plot_signal()
{
  nodes.grainNode.port.onmessage = function(e)
  {
    e.data.buffer1.forEach((value, index, array) =>
    {
      DATA1[index] = Math.abs(value * 128);
    });
    e.data.buffer2.forEach((value, index, array) =>
    {
      DATA2[index] = Math.abs(value * 128);
    });
  }
  
  noFill();  
  strokeWeight(1);
  for(let i = 0; i < SIZE; i++)
  {
    beginShape();
    let posistion = -(frameCount % DEG) * 2 * Math.PI / DEG + i * SPACE;
    let _x = 1 * Math.cos(posistion) * windowWidth/8;
    let _y = 1 * Math.sin(posistion) * windowWidth/8;
    vertex(_x, _y, -DATA1[i]);       
    vertex(_x, _y,  DATA2[i]);    
    if(i == 0)
    {
      stroke(255,0,0);
    }
    else
    {
      stroke(255,255,255);
    }
    endShape(CLOSE);
  } 
}