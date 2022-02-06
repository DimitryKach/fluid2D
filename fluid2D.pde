import java.util.Locale;

// Sim resolution (square)
final int Ndim = 160;
int window_size = 640;
float scale = float(window_size) / float(Ndim);
Fluid fluid;
PImage img;
boolean init = false;
float xoff = 0;
void setup() {
   size(640, 640);
   fluid = new Fluid(0.0001, 0.0005, 0.01);
   //float max = 0;
   //for(int x=0; x<Ndim; x++)
   //{
   //  for(int y=0; y<Ndim; y++)
   //  {
   //    float i = (x-124.5);
   //    float j = (y-124.5);
   //    float dist = sqrt(i*i + j*j);
   //    if(dist > 10) continue;
   //    float dens = (100-dist);
   //    if (dens > max) max=dens;
   //    //if(dist<2) print(dens + "\n");
   //    //print(x + "\n");
   //    //print(y + "\n");
   //    //print(dens + "\n");
   //    fluid.add_dens(int(x/scale)+1,int(y/scale)+1, dens);
   //  }
   //}
   //print(max);
   img = loadImage("D:\\Personal\\processing\\cropped.jpg");
   img.loadPixels();
}

void draw() {
  background(0);
  if(!init)
  {
    //print(1);
    //background(0);
    //for(int x=1; x<=Ndim; x++)
    // {
    //   for(int y=1; y<=Ndim; y++)
    //   {
    //     int p_x = int(x*scale);
    //     int p_y = int(y*scale);
    //     float pixel = red(img.get(p_x, p_y));
    //     //fill(pixel);
    //     fluid.add_dens(x, y, pixel);
    //   }
    // }
    init=true;
  }
  if (mousePressed){
    float dx = (mouseX - pmouseX)*2;
    float dy = (mouseY - pmouseY)*2;
    int i = int(mouseX/scale)+1;
    int j = int(mouseY/scale)+1;
    //print(i, j, "\n");
    for (int x = i-6; x <= i+6; x++){
        for (int y = j-6; y <= j+6; y++){
          if (mouseButton == LEFT) fluid.add_dens(x, y, 45-3*(abs(i-x)+abs(j-y)));
          //if (mouseButton == RIGHT) fluid.add_vel(x, y, dx, dy);
          //print(x, y, fluid.dens[IX(x,y)], " ");
        }
    }
    for (int x = i-3; x <= i+3; x++){
        for (int y = j-3; y <= j+3; y++){
          if (mouseButton == RIGHT) fluid.add_vel(x, y, dx*0.1, dy*0.1);
          //fluid.add_vel(x, y, dx*0.1, dy*0.1);
          //print(x, y, "velocity", dx*0.1, dy*0.1);
        }
    }
  }
  int i = int(320.0/scale)+1;
  int j = int(320.0/scale)+1;
  for (int x = i-5; x <= i+5; x++){
    for (int y = j-5; y <= j+5; y++){
      //fluid.add_dens(x, y, 25);
      //print(x, y, "velocity", dx*0.1, dy*0.1);
    }
  }
  //for (int x = i-5; x <= i+5; x++){
  //  for (int y = j-5; y <= j+5; y++){
  for (int x = 1; x <= Ndim; x++){
    for (int y = 1; y <= Ndim; y++){
      //float dy = map(noise(xoff, 0), 0, 1, -0.75, 0.75);
      float g = 9.8/100;
      //fluid.add_vel(x, y, 0, g*fluid.dens[IX(x, y)]/100);
      //fluid.add_vel(158-x, y-5, -0.75, -dy);
      xoff += 0.5;
      //print(x, y, "velocity", dx*0.1, dy*0.1);
    }
  }
  fluid.step();
  fluid.render_density();
  //fluid.render_velocity();
  noStroke();
  fill(255,200);
  rect(0,0, 60, 20);
  fill(0);
  text("fps: "+String.format(Locale.ENGLISH, "%5.2f", frameRate), 5, 15);
  //saveFrame("capture/cap-######.png");
}
