#include <arduinoFFT.h>

  #include<math.h>
  #define sampling_rate 
  #define pi 3.1415
  
  int LED1=5;
  int LED2=9;
  int LED3=10;
  int LED4=11;
  
  int LED5=6;
  int LED6=2;
  int LED7=3;
  int LED8=4;
  float alpha;
  int t=0;
  int sw;
  float w_1=0.2;
  float w_2=0.3;
  float w_3=0.4;
  float w_4=0.5;
float w_mean=(w_1+w_2+w_3+w_4)/4;
  float alpha_0=0;
  float alpha_1=0;
  float alpha_2=0;
  float alpha_3=0;
  void setup() {
   
    Serial.begin(300);
    pinMode(LED1,OUTPUT);
    pinMode(LED2,OUTPUT);
    pinMode(LED3,OUTPUT);
    pinMode(LED4,OUTPUT);
    pinMode(LED5,OUTPUT);
    pinMode(LED6,OUTPUT);
    pinMode(LED7,OUTPUT);
    pinMode(LED8,OUTPUT);
    ran();
    stable();
     }
  
  int heaviside(double x)
  {
    if(x>=0)
    { 
      return 1;
    }
    else
    {
      return 0;
    }
  }

  void stable()
  {
   while(t<600)
     {
   sw=heaviside(((sin(2*pi*0.35*t))));
   Serial.print(sw);
    digitalWrite(LED1,sw);
    digitalWrite(LED5,sw);
    sw=heaviside(((sin(2*pi*0.35*t))));
    digitalWrite(LED2,sw);
    digitalWrite(LED6,sw);
    sw=heaviside(((sin(2*pi*0.35*t))));
    digitalWrite(LED3,sw);
    digitalWrite(LED7,sw);
    sw=heaviside(((sin(2*pi*0.35*t))));
     digitalWrite(LED4,sw);
     digitalWrite(LED8,sw);
   t=t+2;
  
     }  
   }
  void ran()
  {
  
    while(t<300)
    {
    sw=heaviside(((sin(2*pi*w_1*t))));
    Serial.print(sw);
    digitalWrite(LED1,sw);
    digitalWrite(LED5,sw);
    sw=heaviside((sin(2*pi*w_2*t)));
    digitalWrite(LED2,sw);
    digitalWrite(LED6,sw);
    sw=heaviside((sin(2*pi*w_3*t)));
    digitalWrite(LED3,sw);
    digitalWrite(LED7,sw);
    sw=heaviside((sin(2*pi*w_4*t)));
    digitalWrite(LED4,sw);
    digitalWrite(LED8,sw);
    t=t+1;
    }
  }
  void loop() 
  {
    alpha = pi/2;
 
   
/* 
{
    sw=heaviside(((sin(2*pi*w_1*t+alpha_0))));
    digitalWrite(LED1,sw);
    sw=heaviside((sin(2*pi*w_2*t+alpha_1)));
    digitalWrite(LED2,sw);
    sw=heaviside((sin(2*pi*w_3*t+alpha_2)));
    digitalWrite(LED3,sw);
    sw=heaviside((sin(2*pi*w_4*t+alpha_3)));
    digitalWrite(LED4,sw);
    t=t+100;
 }

 */
 
 //ran();
//stable();
 /*  {
   sw=heaviside(((sin(2*pi*w_mean*t))));
    digitalWrite(LED1,sw);
    sw=heaviside(((sin(2*pi*w_mean*t))));
    digitalWrite(LED2,sw);
    sw=heaviside(((sin(2*pi*w_mean*t))));
    digitalWrite(LED3,sw);
    sw=heaviside(((sin(2*pi*w_mean*t))));
     digitalWrite(LED4,sw);
   t=t+100;
   }
 */  
 
   sw=heaviside((sin(2*pi*0.35*t)));
    digitalWrite(LED5,sw);
   digitalWrite(LED1,sw);
  
   sw=heaviside((sin(2*pi*0.35*t+alpha)));
   digitalWrite(LED6,sw);
   digitalWrite(LED2,sw);
   
   sw=heaviside((sin(2*pi*0.35*t+2*alpha)));
     digitalWrite(LED7,sw);
   digitalWrite(LED3,sw);
 
   sw=heaviside((sin(2*pi*0.35*t+3*alpha)));
   digitalWrite(LED8,sw);
   digitalWrite(LED4,sw);
   
   t=t+300;
  

 
  
   
   
    }

   

  
  
  
