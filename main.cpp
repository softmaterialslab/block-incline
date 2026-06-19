  // This is a code to simulate dynamics of one-particle systems
  // Vikram Jadhao and Vijay G. Jadhao
  // June 19, 2026

  #include<iostream>
  #include<fstream>
  #include<vector>
  #include"vector3d.h"
  #include"particle.h"
  
  const double g = 9.8;
  const int planewidth = 8;
  const int downshift = -4;
  
  using namespace std;
  
  // make movie
  void make_movie(int num, PARTICLE block, PARTICLE& plane, ofstream& outdump)
  {
    
    int L = 24;
    int a, b;
    a = planewidth;
    b = planewidth;
    
    double d = cos(M_PI/4);
    
    outdump << "ITEM: TIMESTEP" << endl;
    outdump << num << endl; // the time at which movie frame is made
    outdump << "ITEM: NUMBER OF ATOMS" << endl;
    outdump << 7 << endl; // hard coded
    outdump << "ITEM: BOX BOUNDS" << endl;
    outdump << -L/2 << "\t" << L/2 << endl;
    outdump << -L/2 << "\t" << L/2 << endl;
    outdump << -L/2 << "\t" << L/2 << endl;
    outdump << "ITEM: ATOMS index type x y" << endl;
    
    PARTICLE P = block;
    PARTICLE Q = block;
    Q.position = block.position + VECTOR3D(d,-d,0);
    PARTICLE R = block;
    R.position = block.position + VECTOR3D(d,d,0);
    PARTICLE S = block;
    S.position = block.position + VECTOR3D(2*d,0,0);
    
    //outdump << 1 << "   " << block.ty << "   " << block.position.x << "   " << block.position.y << endl;
    
    outdump << 1 << "   " << P.ty << "   " << P.position.x << "   " << P.position.y << endl;
    outdump << 2 << "   " << Q.ty << "   " << Q.position.x << "   " << Q.position.y << endl;
    outdump << 3 << "   " << R.ty << "   " << R.position.x << "   " << R.position.y << endl;
    outdump << 4 << "   " << S.ty << "   " << S.position.x << "   " << S.position.y << endl;
    
    PARTICLE A = plane;
    PARTICLE B = plane;
    B.position = plane.position + VECTOR3D(0,b,0);
    PARTICLE C = plane;
    C.position = plane.position + VECTOR3D(a,0,0);
    
    //outdump << 2 << "   " << plane.ty << "   " << plane.position.x << "   " << plane.position.y << endl;
    
    outdump << 5 << "   " << A.ty << "   " << A.position.x << "   " << A.position.y << endl;
    outdump << 6 << "   " << B.ty << "   " << B.position.x << "   " << B.position.y << endl;
    outdump << 7 << "   " << C.ty << "   " << C.position.x << "   " << C.position.y << endl;
  }

  void compute_energy(PARTICLE& block, PARTICLE& plane, double theta)
  {
      // pe of the block
      block.pe = block.m * g * block.position.y;
      plane.pe = 0;
  }

  void compute_force(PARTICLE& block, PARTICLE& plane, double theta)
  {
      double r = block.m/plane.m;   // ratio of block mass and plane mass
      r = 0;
          
      // force for the block 
      block.force.x = block.m * g * cos(theta) * sin(theta) / (1 + r * sin(theta) * sin(theta));
      block.force.y = -block.m * g * (1 + r) * sin(theta) * sin(theta) / (1 + r * sin(theta) * sin(theta));
      
     
      // force on the plane
      plane.force.x = -plane.m * r * g * cos(theta) * sin(theta) / (1 + r * sin(theta) * sin(theta)); 
  }

  int main(int argc, char* argv[])
  {
      cout << "Simulating block plane dynamics" << endl;

      cout << "Set up initial conditions" << endl;

      double m = 1;
      double M = 1;

      //double k = 1; // spring constant

      //cout << "enter theta" << endl;
      double theta = M_PI / 4;
      //cin >> theta;

      // simulating the block-spring system

      PARTICLE block = PARTICLE(m,1,VECTOR3D(1,planewidth-1+downshift,0),VECTOR3D(0,0,0));
      
      PARTICLE plane = PARTICLE(M,2,VECTOR3D(0,downshift,0),VECTOR3D(0,0,0));
      
      // initial force on the block
      compute_force(block,plane,theta);

      cout << "mass of the block is " << block.m << endl;
      cout << "initial position of the block is " << block.position.x << endl;
      cout << "initial velocity of the block is " << block.velocity.x << endl;
      cout << "initial force on the block is " << block.force.x << endl;
      cout << "inclined plane angle " << theta * 180 / M_PI << endl;

      block.kinetic_energy();
      plane.kinetic_energy();
      compute_energy(block,plane,theta);

      // writing files to store data
      char filename_block[200];
      sprintf(filename_block, "energyMinf.out");
      ofstream block_result(filename_block, ios::out);
      
      ofstream list_particles("movieMinf.out", ios::out); // create a file to store and visualize 3D data

      // start by recording the initial condition (state) of the system
      block_result << 0 << "  " << block.ke << "  " << block.pe << "  " << block.ke + block.pe << "  " << plane.ke << "  " << plane.pe << "  " << plane.ke + plane.pe << "  " << block.ke + block.pe + plane.ke + plane.pe << endl;

      make_movie(0,block,plane,list_particles);
      
      double t = 2; 	// total duration of time for which dynamics is desired
      int S = 2000;
      double dt = t/S; // dt = 20/20000 = 0.001
      
      int freq = 10;

      // simulation
      for (int num = 1; num < S; num++)
      {
          // update velocity (half step)
          block.update_velocity(dt/2);
          plane.update_velocity(dt/2);

          // update position
          block.update_position(dt);
          plane.update_position(dt);

          // recalculate force
          compute_force(block,plane,theta);

          // update velocity again
          block.update_velocity(dt/2);
          plane.update_velocity(dt/2);

          // useful to collect energy information
          block.kinetic_energy();
          plane.kinetic_energy();
          compute_energy(block,plane,theta);

          // file the simulation updates
          block_result << num*dt << "  " << block.ke << "  " << block.pe << "  " << block.ke + block.pe << "  " << plane.ke << "  " << plane.pe << "  " << plane.ke + plane.pe << "  " << block.ke + block.pe + plane.ke + plane.pe << endl;

          if (num%freq == 0)
            make_movie(num,block,plane,list_particles);
      }

      return 0;
  }
