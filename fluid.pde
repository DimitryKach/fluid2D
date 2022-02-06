int IX(int i, int j){
  i = constrain(i, 0, Ndim+1);
  j = constrain(j, 0, Ndim+1);
  return i + j * (Ndim+2);
}

void SWAP(float[] x0, float[] x){
    float[] temp=new float[x0.length];
    arrayCopy(x0, temp);
    arrayCopy(x, x0);
    arrayCopy(temp, x);
}

class Fluid{
   // Following the Jos Stam paper, we allocate arrays (single dimension) for velocity and for density
  int size;
  float diff;
  float visc;
  float step_size;
  float[] u;
  float[] v;
  float[] u_prev;
  float[] v_prev;
  float[] dens;
  float[] dens_prev;
  // Sources array
  float[] s;
  
  Fluid(float diff, float visc, float step_size) {
    this.size = (Ndim+2) * (Ndim+2);
    this.step_size = step_size;
    this.diff = diff;
    this.visc = visc;
    this.u = new float[this.size];
    this.v = new float[this.size];
    this.u_prev = new float[this.size];
    this.v_prev = new float[this.size];
    this.dens = new float[this.size];
    this.dens_prev = new float[this.size];
    // Sources array
    this.s = new float[this.size];
  }
  
  void add_dens( int i, int j, float d_dens){
      this.dens[IX(i, j)] += d_dens;
  }
  
  void add_vel( int i, int j, float d_vel_x, float d_vel_y){
      this.u[IX(i,j)] += d_vel_x;
      this.v[IX(i,j)] += d_vel_y;
  }
  
  void render_density(){
      for (int i=1; i<=Ndim; i++){
          for (int j=1; j<=Ndim; j++){
              float x = i*scale;
              float y = j*scale;
              noStroke();
              PVector col = new PVector(0,0,0);
              float dens = this.dens[IX(i,j)];
              float max = 1000;
              //if(this.dens[IX(i,j)]!=0){
              //  dens = min(max, this.dens[IX(i,j)]);
              //}
              rainbow(dens/max, col);
              fill(255, this.dens[IX(i,j)]);
              //fill(col.x, col.y, col.z, 255);
              //fill(255, 0);
              rect(x, y, scale, scale);
          }
      }
  }
  
  void rainbow(float d, PVector RGB)
  {
     float dx=1;
     float value = d;
     value = value/0.2;
     if (value<0) value=0; if (value>1) value=1;
     if (value==0) return;
     value = (6-2*dx)*d+dx;

     RGB.x = max(0.0,(3-abs(value-4)-abs(value-5))/2)*255;
     RGB.y = max(0.0,(4-abs(value-2)-abs(value-4))/2)*255;
     RGB.z = max(0.0,(3-abs(value-1)-abs(value-2))/2)*255; 
  }
  
  void render_velocity(){
    for (int j = 1; j <= Ndim; j++) {
      for (int i = 1; i <= Ndim; i++) {
        float x = i * scale;
        float y = j * scale;
        stroke(255);
        strokeWeight(1);
        PVector v = new PVector(this.u[IX(i, j)], this.v[IX(i, j)]);
        //print(v);
        if (v.mag() > 0.1){
          float max = 1000;
          if (v.mag() > max){
            v.normalize();
            v = v.mult(max);
          }
          line(x, y, x + v.x, y + v.y);
        }
      }
    }
  }
  
  void _diffuse(float[] x, float[] x0, float diff, float dt, int b){
    // we are trying to diffuse our density into nearby cells. Going straight ahead from the current place is possible, but is unstable (Forward Euler).
    // But lets still write it out. If we imagine just diffusing in 1 direction (hence, only 1 dimension), we can think of maybe just 5 cells that we want to work with
    // we consider the indices from 0-4. Indices 0 and 4 are boundaries, so we don't care to diffuse there. So lets then name cell ids as x_i, i.e. x_0, x_1, ... x_n
    // and x0_0, x0_1, ... x0_n.
    // We need to remember that we also have nearby cells, and they also diffuse into the current cell. And our current cell will also loose density to nearby cells.
    // So lets describe how the next value of, say cell x_2 would be computed. We first will loose half of x_2's density, so x_2 = x0_2 - a*(2*x0_2). And then it will have
    // the nearby densities added: x_2 = x0_2 + a*(x0_1 + x0_3) - a*(2*x0_2). Lets simplify: x_2 = x0_2 + a*(x0_1 + x0_3 - 2*x0_2)
    // But lets make that into a matrix form. For that, lets expand the equation, and isolate the variables best we can.
    // x_2 = x0_2 + a*x0_1 + a*x0_3 - 2*a*x0_2
    // x_2 = x0_1*a + x0_2*(1-2*a) + x0_3*a
    // If we look at the last equation we got, we have a nice looking trio of constant that multiply against each variable, and then they are all added up. This looks like
    // the good old dot product of 2 vectors! Lets put all our contants into a vector we call A = {a, 1-2*a, a}, the current vector X0s(for short, since we are only using
    // a part of our full vector) will be X0s = {x0_1, x0_2, x0_3}. So x_2 = dot(A, X0s)! Well that process can be also put into a matrix, so that we can cover our entire
    // vector X0 with size 5. First, in vector form, we can generalize, that, x_i = dot(A, {x0_(i-1), x0_i, x0_(i+1)}).
    // In matrix form, ignoring index 0, we form an A matrix that looks like this:
    // 0   0   0   0   0
    // a  1-2a a   0   0
    // 0   a  1-2a a   0
    // 0   0   a  1-2a a
    // 0   0   0   0   0
    // And if now think of the whole process, our new vector x = Ax0.
    // But our story doesn't end here... Since our A depends on dt, lets call the current matrix A, A(dt), since we are going forward in time step.
    // This process is just basic forward euler integration. However, it is known to be unstable, so what Jos Stam proposes, is that we use backward euler. Here is how
    // we can look at it. If our x was at t, and x0 is same as x at t-1, we can re-write the formula as x(t) = A(dt)x(t-1). But lets imagine that we know our current
    // density, and we would rather find out what was the previous density that made this current one what it is. So lets flip some signs - x(t-1) = A(-dt)x(t). So
    // what happened here is that since we *know* our current density, we need to step BACK in time (hence the -dt) and run our euler backwards. To understand what A(-dt)
    // would look like, we need to quickle expand this around a single variable again. Lets take x0_2, for example - x0_2 = x_2 - a*(x_1 + x_3 - 2*x_2).
    // Since we are going backwards, instead of adding our estimation, I need to remove it, hence the main "-" operation. Going back to the expanded dot product version
    // we get x0_2 = x_2(1+2a)-x_3a-x_1*a, where vector a is {-a, 1+2a, -a}. If we look back at our matrix, and replace the previous values with these, we will get our
    // A(-dt). And here comes the cool part. Since in reality we actually know our PREVIOUS densities at x(t-1), and not at x(t), the formula A(-dt)x(t)=x(t-1) takes
    // the shape of a classical Ax=b linear system! All we need to do is solve for x. In Jos Stems paper he is solving that using the Gauss Seidel relaxation.
    // Its Gauss Seidel because we are succesively updating our values as we find them. Lets derive the formula below:
    // We have a varibale "a" which is our diffuse rate. Given our 2D grid of cells, we need to remember we diffuse each cell into 4 other cells. So lets rewrite our
    // previous equation here but with 4 elements - two on the X axis, and two on the Y (given by i and j):
    // x0_[i,j] = x_[i,j] - a*(x_[i-1,j] + x_[i+1,j] + x_[i,j-1] + x_[i,j+1] - 4*x_[i,j])
    // x0_[i,j] = x_[i,j](1+4a) - ax_[i-1,j] - ax_[i+1,j] - ax_[i,j-1] - ax_[i,j+1]
    // Now lets solve for x_[i,j]:
    // x_[i,j](1+4a) = x0_[i,j] + ax_[i-1,j] + ax_[i+1,j] + ax_[i,j-1] + ax_[i,j+1]
    // x_[i,j] = (x0_[i,j] + ax_[i-1,j] + ax_[i+1,j] + ax_[i,j-1] + ax_[i,j+1])/(1+4a)
    // x_[i,j] = (x0_[i,j] + a(x_[i-1,j] + x_[i+1,j] + x_[i,j-1] + x_[i,j+1]))/(1+4a) -> our final iterator function!
    // And, as we iterate, we overwrite our previous guesses in x, and thus create what is known as Gauss-Seidel numerical method.
    float a = dt * diff * Ndim * Ndim;
    float c = (1+4*a);
    boolean w_sor = false; 
    int iter = 20;
    if (w_sor){
      for (int k=0; k<iter; k++){
        for (int i=1; i<=Ndim; i++){
          for (int j=1; j<=Ndim; j++){
            float w = 1.3;
            x[IX(i,j)] = (1-w)*x[IX(i,j)] + w*((x0[IX(i,j)] + a * (x[IX(i-1,j)] + x[IX(i+1,j)] + x[IX(i,j-1)] + x[IX(i,j+1)])) / (1+4*a));
          }
        }
      }
    }
    else{
      gauss_seidel(Ndim, x, x0, a, c, iter, b);
    }
  }
  
  void _advect(float[] x, float[] x0, float[] u, float[] v, float dt){
    float dt0 = dt * Ndim;
    for (int i=1; i<=Ndim; i++){
      for (int j=1; j<=Ndim; j++){
        float vx = u[IX(i,j)];
        float vy = v[IX(i,j)];
        //print(vx, vy, "\n");
        // We need to interpolate, but, we need to do so from the center
        //float x = i + 0.5;
        //float y = j + 0.5;
        // step backwards
        float xi = i - dt0*vx;
        float yj = j - dt0*vy;
        // Lets make sure we stay within the bounds...
        xi = max(0.5, min(float(Ndim)+0.5, xi));
        yj = max(0.5, min(float(Ndim)+0.5, yj));
        // Lets determin the new cell
        int i0 = int(xi);
        int j0 = int(yj);
        int i1 = i0+1;
        int j1 = j0+1;
        // bilinear interpolation
        float q11 = x0[IX(i0, j0)];
        float q12 = x0[IX(i1, j0)];
        float q21 = x0[IX(i0, j1)];
        float q22 = x0[IX(i1, j1)];
        
        float x_interp0 = xi - i0;
        float x_interp1 = 1.0 - x_interp0;
        float y_interp0 = yj - j0;
        float y_interp1 = 1.0 - y_interp0;
        
        float fx0_y0 = q11*x_interp1 + q12*x_interp0;
        float fx0_y1 = q21*x_interp1 + q22*x_interp0;
        
        float fxy = fx0_y0*y_interp1 + fx0_y1*y_interp0;
        
        x[IX(i,j)] = fxy;
      }
    }
  }
  
  void dens_step(){
    SWAP(this.dens_prev, this.dens);
    this._diffuse(this.dens, this.dens_prev, this.diff, this.step_size, 0);
    SWAP(this.dens_prev, this.dens);
    this._advect(this.dens, this.dens_prev, this.u, this.v, this.step_size);
  }
  
  void vel_step(){
    SWAP(this.u_prev, this.u); SWAP(this.v_prev, this.v);
    this._project(this.u, this.v, this.u_prev, this.v_prev);
    SWAP(this.u_prev, this.u); SWAP(this.v_prev, this.v);
    this._diffuse(this.u, this.u_prev, this.visc, this.step_size, 1);
    this._diffuse(this.v, this.v_prev, this.visc, this.step_size, 2);
    SWAP(this.u_prev, this.u); SWAP(this.v_prev, this.v);
    this._project(this.u, this.v, this.u_prev, this.v_prev);
    SWAP(this.u_prev, this.u); SWAP(this.v_prev, this.v);
    this._advect(this.u, this.u_prev, this.u_prev, this.v_prev, this.step_size);
    this._advect(this.v, this.v_prev, this.u_prev, this.v_prev, this.step_size);
    SWAP(this.u_prev, this.u); SWAP(this.v_prev, this.v);
    this._project(this.u, this.v, this.u_prev, this.v_prev);
  }
  
  boolean is_solid(int i, int j)
  {
    return false;
  }
  
  void _project(float[] u, float[] v, float[] u0, float[] v0){
    // Here we need to solve for incompressibility of the fluid. As described in the Jos Stam paper, we are solving the poisson equation. We arrive at that equation
    // by using the result of the Hodge decompositions (that portion I don't know how that works). The main result of that decomposition states that a velocity field
    // is the sum of a mass conserving field and a gradient field. So if we write that out, V = V' + G(f), where V is the velocity field, V' is the mass conserving
    // field, and G(f) is the gradient field. We are after the V', or the mass conserving field. That makes our function look like V' = V - G(f).
    // In our case, the field we are taking the gradient of is the pressure field. This makes our function look like V' = V - G(p).
    // Being that we know that V' is incompressible, but also, mass conserving, which means:
    // D(V') = 0
    // which states that the divergence of the V' field is equal to 0. Using that information, lets look at what we can derive:
    // V' = V-G(p)
    // DV' = D(V-G(p)) = 0
    // DV - DG(p) = 0
    // DV = DG(p)
    // What we now have is an interesting form. Lets remember that the Divergence of a Gradient is actually a Laplacian, so lets call DG = L.
    // The Poisson equation states L(phi) = f, where phi and f are real or compelx valued fanctions (on a manifold, although this is where I get lost).
    // If we look at our equation, and realize that since V is a vector field, DV makes is a scalar field. P is also a scalar field. So lets call DV=f, we get
    // L(p)=f - our Poisson equation! So all we need to do now is to solve it, and Jos Stam uses the Gauss Seidel once again, with his pressures having an initial
    // guess of 0. This will give us the p. Once we have the p, we can go back to the original field equation, and solve for V'!
    // So lets outline our steps:
    // 1) For every grid cell (from 1-N, leaving the border grids out of this) we want the divergence of each cell, and then set that cells pressure to 0
    float[] p = new float[u.length];
    float[] div = new float[u.length];
    float h = 1.0/Ndim;
    for (int i=1; i<=Ndim; i++){
      for (int j=1; j<=Ndim; j++){
        // The divergence here is a sum of in and out flow of fluid into the cell. We need to set it to negative to account for the in flow (since divergence is
        // out by default, so positivity is out-flow, while negativity is in-flow)
        div[IX(i,j)] = -0.5*h*(u0[IX(i+1,j)] - u0[IX(i-1,j)] + v0[IX(i,j+1)] - v0[IX(i,j-1)]);
        p[IX(i,j)] = 0.0;
      }
    }
    set_bnd(Ndim, 0, div); set_bnd(Ndim, 0, p);
    // 2) Lets use Gauss-Seidel to approximate the value of p. We do that by iteratively adding the divergence at the current cell and the pressures of the
    //    surrounding cells and divided by 4 (since I guess there are 4 neighbors...
    //    So, p_xy = (div_xy + p_x-1y + p_x+1y + p_xy-1 + p_xy+1)/4;
    int iter = 20;
    float a = 1.0;
    float c = 4.0;
    gauss_seidel(Ndim, p, div, a, c, iter, 0);
    arrayCopy(u0, u);
    arrayCopy(v0, v);
    // 3) Finally, lets update our velocity
    for (int i=1; i<=Ndim; i++){
      for (int j=1; j<=Ndim; j++){
        u[IX(i,j)] -= 0.5*(p[IX(i+1,j)] - p[IX(i-1,j)])/h;
        v[IX(i,j)] -= 0.5*(p[IX(i,j+1)] - p[IX(i,j-1)])/h;
      }
    }
    set_bnd(Ndim, 1, u); set_bnd(Ndim, 2, v);
  };
  
  void step(){
    this.vel_step();
    this.dens_step();
  }
}

void gauss_seidel( int N, float[] x, float[] x0, float a, float c, int iter, int b){
  for (int k=0; k<iter; k++){
      for (int i=1; i<=N; i++){
        for (int j=1; j<=N; j++){
          x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)] + x[IX(i+1,j)] + x[IX(i,j-1)] + x[IX(i,j+1)])) / c;
        }
      }
      set_bnd(N, b, x);
    }
}

void set_bnd ( int N, int b, float[] x )
{
  int i;
  for ( i=1 ; i<=N ; i++ ) {
    x[IX(0 ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
    x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
    x[IX(i,0 )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
    x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
  }
  x[IX(0 ,0 )] = 0.5*(x[IX(1,0 )]+x[IX(0 ,1)]);
  x[IX(0 ,N+1)] = 0.5*(x[IX(1,N+1)]+x[IX(0 ,N )]);
  x[IX(N+1,0 )] = 0.5*(x[IX(N,0 )]+x[IX(N+1,1)]);
  x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N )]);
}
