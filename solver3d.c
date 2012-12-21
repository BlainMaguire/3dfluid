 /* 

A 3D Real Time Fluid Solver based on Jos Stam's fluid solver (stable N-S solver).
 
Reference: Jos Stam, "Real-Time Fluid Dynamics for Games". Proceedings of the Game Developer Conference, March 2003.
  *http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf]

*/

#define IX(i,j,k) ((i)+(M+2)*(j) + (M+2)*(N+2)*(k))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define MAX(a,b)            (((a) > (b)) ? (a) : (b))
#define LINEARSOLVERTIMES 10

void add_source ( int M, int N, int O, float * x, float * s, float dt )
{
	int i, size=(M+2)*(N+2)*(O+2);
	for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}

void  set_bnd ( int M, int N, int O, int b, float * x )
{

        // bounds are cells at faces of the cube

	int i, j;

	for ( i=1 ; i<=M ; i++ ) {
		for ( j=1 ; j<=N ; j++ ) {
			x[IX(i,j,0 )] = b==3 ? -x[IX(i,j,1)] : x[IX(i,j,1)];
			x[IX(i,j,O+1)] = b==3 ? -x[IX(i,j,O)] : x[IX(i,j,O)];
		}
	}

    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=O ; j++ ) {
            x[IX(0  ,i, j)] = b==1 ? -x[IX(1,i,j)] : x[IX(1,i,j)];
            x[IX(M+1,i, j)] = b==1 ? -x[IX(M,i,j)] : x[IX(M,i,j)];
        }
    }

    for ( i=1 ; i<=M ; i++ ) {
        for ( j=1 ; j<=O ; j++ ) {
            x[IX(i,0,j )] = b==2 ? -x[IX(i,1,j)] : x[IX(i,1,j)];
            x[IX(i,N+1,j)] = b==2 ? -x[IX(i,N,j)] : x[IX(i,N,j)];
        }
    }

    x[IX(0  ,0, 0  )] = 1.0/3.0*(x[IX(1,0,0  )]+x[IX(0  ,1,0)]+x[IX(0 ,0,1)]);
    x[IX(0  ,N+1, 0)] = 1.0/3.0*(x[IX(1,N+1, 0)]+x[IX(0  ,N, 0)] + x[IX(0  ,N+1, 1)]);

    x[IX(M+1,0, 0 )] = 1.0/3.0*(x[IX(M,0,0  )]+x[IX(M+1,1,0)] + x[IX(M+1,0,1)]) ;
    x[IX(M+1,N+1,0)] = 1.0/3.0*(x[IX(M,N+1,0)]+x[IX(M+1,N,0)]+x[IX(M+1,N+1,1)]);

    x[IX(0  ,0, O+1 )] = 1.0/3.0*(x[IX(1,0,O+1  )]+x[IX(0  ,1,O+1)]+x[IX(0 ,0,O)]);
    x[IX(0  ,N+1, O+1)] = 1.0/3.0*(x[IX(1,N+1, O+1)]+x[IX(0  ,N, O+1)] + x[IX(0  ,N+1, O)]);

    x[IX(M+1,0, O+1 )] = 1.0/3.0*(x[IX(M,0,O+1  )]+x[IX(M+1,1,O+1)] + x[IX(M+1,0,O)]) ;
    x[IX(M+1,N+1,O+1)] = 1.0/3.0*(x[IX(M,N+1,O+1)]+x[IX(M+1,N,O+1)]+x[IX(M+1,N+1,O)]);
}

void lin_solve ( int M, int N, int O, int b, float * x, float * x0, float a, float c )
{
	int i, j, k, l;
        
	// iterate the solver
	for ( l=0 ; l<LINEARSOLVERTIMES ; l++ ) {
		// update for each cell
		for ( i=1 ; i<=M ; i++ ) { for ( j=1 ; j<=N ; j++ ) { for ( k=1 ; k<=O ; k++ ) {
			x[IX(i,j,k)] = (x0[IX(i,j,k)] + a*(x[IX(i-1,j,k)]+x[IX(i+1,j,k)]+x[IX(i,j-1,k)]+x[IX(i,j+1,k)]+x[IX(i,j,k-1)]+x[IX(i,j,k+1)]))/c;
		}}}
        set_bnd ( M, N, O, b, x );
	}
}

void diffuse (  int M, int N, int O, int b, float * x, float * x0, float diff, float dt )
{
	int max = MAX(MAX(M, N), MAX(N, O));
	float a=dt*diff*max*max*max;
	lin_solve ( M, N, O, b, x, x0, a, 1+6*a );
}

void advect (  int M, int N, int O, int b, float * d, float * d0, float * u, float * v, float * w, float dt )
{
	int i, j, k, i0, j0, k0, i1, j1, k1;
	float x, y, z, s0, t0, s1, t1, u1, u0, dtx,dty,dtz;
	
	dtx=dty=dtz=dt*MAX(MAX(M, N), MAX(N, O));

	for ( i=1 ; i<=M ; i++ ) { for ( j=1 ; j<=N ; j++ ) { for ( k=1 ; k<=O ; k++ ) {
		x = i-dtx*u[IX(i,j,k)]; y = j-dty*v[IX(i,j,k)]; z = k-dtz*w[IX(i,j,k)];
		if (x<0.5f) x=0.5f; if (x>M+0.5f) x=M+0.5f; i0=(int)x; i1=i0+1;
		if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
		if (z<0.5f) z=0.5f; if (z>O+0.5f) z=O+0.5f; k0=(int)z; k1=k0+1;

		s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1; u1 = z-k0; u0 = 1-u1;
		d[IX(i,j,k)] = s0*(t0*u0*d0[IX(i0,j0,k0)]+t1*u0*d0[IX(i0,j1,k0)]+t0*u1*d0[IX(i0,j0,k1)]+t1*u1*d0[IX(i0,j1,k1)])+
			s1*(t0*u0*d0[IX(i1,j0,k0)]+t1*u0*d0[IX(i1,j1,k0)]+t0*u1*d0[IX(i1,j0,k1)]+t1*u1*d0[IX(i1,j1,k1)]);
	}}}
    
    set_bnd (M, N, O, b, d );
}

void project ( int M, int N, int O, float * u, float * v, float * w, float * p, float * div )
{
	int i, j, k;

	for ( i=1 ; i<=M ; i++ ) { for ( j=1 ; j<=N ; j++ ) { for ( k=1 ; k<=O ; k++ ) {
		div[IX(i,j,k)] = -1.0/3.0*((u[IX(i+1,j,k)]-u[IX(i-1,j,k)])/M+(v[IX(i,j+1,k)]-v[IX(i,j-1,k)])/M+(w[IX(i,j,k+1)]-w[IX(i,j,k-1)])/M);
		p[IX(i,j,k)] = 0;
	}}}	
	
	set_bnd ( M, N, O, 0, div ); set_bnd (M, N, O, 0, p );

	lin_solve ( M, N, O, 0, p, div, 1, 6 );

	for ( i=1 ; i<=M ; i++ ) { for ( j=1 ; j<=N ; j++ ) { for ( k=1 ; k<=O ; k++ ) {
		u[IX(i,j,k)] -= 0.5f*M*(p[IX(i+1,j,k)]-p[IX(i-1,j,k)]);
		v[IX(i,j,k)] -= 0.5f*M*(p[IX(i,j+1,k)]-p[IX(i,j-1,k)]);
		w[IX(i,j,k)] -= 0.5f*M*(p[IX(i,j,k+1)]-p[IX(i,j,k-1)]);
	}}}
	
	set_bnd (  M, N, O, 1, u ); set_bnd (  M, N, O, 2, v );set_bnd (  M, N, O, 3, w);
}

void dens_step ( int M, int N, int O, float * x, float * x0, float * u, float * v, float * w, float diff, float dt )
{
	add_source ( M, N, O, x, x0, dt );
	SWAP ( x0, x ); diffuse ( M, N, O, 0, x, x0, diff, dt );
	SWAP ( x0, x ); advect ( M, N, O, 0, x, x0, u, v, w, dt );
}

void vel_step ( int M, int N, int O, float * u, float * v,  float * w, float * u0, float * v0, float * w0, float visc, float dt )
{
	add_source ( M, N, O, u, u0, dt ); add_source ( M, N, O, v, v0, dt );add_source ( M, N, O, w, w0, dt );
	SWAP ( u0, u ); diffuse ( M, N, O, 1, u, u0, visc, dt );
	SWAP ( v0, v ); diffuse ( M, N, O, 2, v, v0, visc, dt );
	SWAP ( w0, w ); diffuse ( M, N, O, 3, w, w0, visc, dt );
	project ( M, N, O, u, v, w, u0, v0 );
	SWAP ( u0, u ); SWAP ( v0, v );SWAP ( w0, w );
	advect ( M, N, O, 1, u, u0, u0, v0, w0, dt ); advect ( M, N, O, 2, v, v0, u0, v0, w0, dt );advect ( M, N, O, 3, w, w0, u0, v0, w0, dt );
	project ( M, N, O, u, v, w, u0, v0 );
}
