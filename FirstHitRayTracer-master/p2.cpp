#include <GL/glut.h>
#include<iostream>
#include<math.h>

struct ray
{ float x,y,z;
};

struct rgb_values
{
    float r; float g; float b;
};

float pos=1.0;
float delta=0.0001;

void display()
{
    int i,j,n_spec_power;
    //float infinity;
    float disc,beta,gamma,max,costheta1,costheta2;
    float ddotesubc, ddotd, e_cdote_c_minus_rsquare, rsquare,ndotl,ndoth, NdotLdiffuse;
    ray a,b,c,d,e,ct1,ct2,ia,i_ds,l,n,pt1,pt2,h,i_sp, A1, B1, C1, D1, N1, N2,E1,N ;
    float A,B,C,D,E,F,G,H,I,J,K,L,M,ka,kd,t1,t2,deno,ks,min,t3,t4,t5,t6,t7,t8;
    rgb_values pixels[512][512],l_ambient, l_ds, l_shading, l_sp,L_diffuse,v1,v2,v3;
    
    
    
    
    
    //specular shading
    ks=0.8;
    n_spec_power=1.1;
    
    i_sp.x=1.5;
    i_sp.y=1.5;
    i_sp.z=1.5;
    
    
    //diffuse shading
    kd=0.9;
    
    i_ds.x=0.8;
    i_ds.y=0.8;
    i_ds.z=0.8;
    
    
    
    
    
    l.x=0.0;
    l.y=-1.0;
    l.z=0.0;
    
    
    
    //ambient shading
    
    ka=0.4;
    ia.x=0.5;
    ia.y=0.5;
    ia.z=0.5;
    
    l_ambient.r=ka*ia.x;
    l_ambient.g=ka*ia.y;
    l_ambient.b=ka*ia.z;
    
    
    
    glClearColor(0.0, 0.0, 0.0, 1.0);
    
    
    
    d.x=pos;
    d.y=0;
    d.z=0;
    
    
    
    
    
    
    e.x=0.0;
    
    a.x=150.0;
    a.y=150.0;
    a.z=200.0;
    
    b.x=150.0;
    b.y=275.0;
    b.z=300.0;
    
    ct1.x=150.0;
    ct1.y=275.0;
    ct1.z=150.0;
    
    ct2.x=-50.0;
    ct2.y=175.0;
    ct2.z=200.0;
    
    min=40000;
    v1.r=1.0;
    v1.g=0.0;
    v1.b=0.0;
    
    v2.r=0.0;
    v2.g=0.0;
    v2.b=1.0;
    
    v3.r=1.0;
    v3.g=1.0;
    v3.b=0.0;
    
    
    for(i=0;i<512;i++)
    {
        for(j=0;j<512;j++)
        {
            e.y=i;
            e.z=j;
            
            rsquare=2500;
            
            c.x=90.0;
            c.y=90.0;
            c.z=90.0;
            
            ray e_sub_c;
            e_sub_c.x=e.x-c.x;
            e_sub_c.y=e.y-c.y;
            e_sub_c.z=e.z-c.z;
            
            ray d_dot_esubc;
            d_dot_esubc.x=d.x*e_sub_c.x;
            d_dot_esubc.y=d.y*e_sub_c.y;
            d_dot_esubc.z=d.z*e_sub_c.z;
            
            ddotesubc=d_dot_esubc.x+d_dot_esubc.y+d_dot_esubc.z;
            
            
            ray d_dot_d;
            d_dot_d.x=d.x*d.x;
            d_dot_d.y=d.x*d.y;
            d_dot_d.z=d.x*d.z;
            
            ddotd=d_dot_d.x+d_dot_d.y+d_dot_d.z;
            
            ray e_sub_cdote_sub_c;
            e_sub_cdote_sub_c.x=e_sub_c.x*e_sub_c.x;
            e_sub_cdote_sub_c.y=e_sub_c.y*e_sub_c.y;
            e_sub_cdote_sub_c.z=e_sub_c.z*e_sub_c.z;
            
            e_cdote_c_minus_rsquare=e_sub_cdote_sub_c.x+e_sub_cdote_sub_c.y+e_sub_cdote_sub_c.z-rsquare;
            
            
            disc=(ddotesubc*ddotesubc)-(ddotd*e_cdote_c_minus_rsquare);
            if(disc<0)
            {
                
            }
            else
            {
                t1=(-(ddotesubc)+sqrt(disc))/ddotd;
                t2=(-(ddotesubc)-sqrt(disc))/ddotd;
                
                
                
                
                if(t1>t2)
                {
                    min=t2;
                    pt1.x=e.x+(t2*d.x);
                    pt1.y=e.y+(t2*d.y);
                    pt1.z=e.z+(t2*d.z);
                    
                }
                else
                {
                    min=t1;
                    pt1.x=e.x+(t1*d.x);
                    pt1.y=e.y+(t1*d.y);
                    pt1.z=e.z+(t1*d.z);
                    
                }
                n.x=pt1.x-c.x;
                n.y=pt1.y-c.y;
                n.z=pt1.z-c.z;
                
                deno=sqrt((n.x*n.x)+(n.y*n.y)+(n.z*n.z));
                n.x=n.x/deno;
                n.y=n.y/deno;
                n.z=n.z/deno;
                
                ray n_dot_l;
                n_dot_l.x=n.x*l.x;
                n_dot_l.y=n.y*l.y;
                n_dot_l.z=n.z*l.z;
                
                ndotl=n_dot_l.x+n_dot_l.y+n_dot_l.z;
                
                if(ndotl>0.0)
                {
                    max=ndotl;
                }
                else
                {
                    max=0.0;
                }
                
                l_ds.r=kd*i_ds.x*max;
                l_ds.g=kd*i_ds.y*max;
                l_ds.b=kd*i_ds.z*max;
                
                
                h.x=d.x+l.x;
                h.y=d.y+l.y;
                h.z=d.z+l.z;
                
                deno=sqrt((h.x*h.x)+(h.y*h.y)+(h.z*h.z));
                h.x=h.x/deno;
                h.y=h.y/deno;
                h.z=h.z/deno;
                
                
                
                
                
                ray n_dot_h;
                n_dot_h.x=n.x*h.x;
                n_dot_h.y=n.y*h.y;
                n_dot_h.z=n.z*h.z;
                
                ndoth=n_dot_h.x+n_dot_h.y+n_dot_h.z;
                
                if(ndoth>0.0)
                {
                    max=ndoth;
                }
                else
                {
                    max=0.0;
                }
                
                l_sp.r=ks*i_sp.x*pow(max,n_spec_power);
                l_sp.g=ks*i_sp.y*pow(max,n_spec_power);
                l_sp.b=ks*i_sp.z*pow(max,n_spec_power);
                
                
                
                
                
                
                
                
                
                /*l_shading.r=l_sp.r;
                 l_shading.g=l_sp.g;
                 l_shading.b=l_sp.b;*/
                l_shading.r=l_ds.r+l_ambient.r+l_sp.r;
                l_shading.g=l_ds.g+l_ambient.g+l_sp.g;
                l_shading.b=l_ds.b+l_ambient.b+l_sp.b;
                
                
                /*pixels[i][j].r=1.0*l_ambient.r;
                 pixels[i][j].g=0.0*l_ambient.g;
                 pixels[i][j].b=0.0*l_ambient.b;
                 
                 
                 pixels[i][j].r=1.0;
                 pixels[i][j].g=0.0;
                 pixels[i][j].b=0.0; */
                
                pixels[i][j].r=1.0*l_shading.r;
                pixels[i][j].g=0.0*l_shading.g;
                pixels[i][j].b=0.0*l_shading.b;
            }
            
            
            
            rsquare=3600;
            
            e.y=i;
            e.z=j;
            c.x=400.0;
            c.y=400.0;
            c.z=400.0;
            
            
            
            e_sub_c.x=e.x-c.x;
            e_sub_c.y=e.y-c.y;
            e_sub_c.z=e.z-c.z;
            
            //    ray d_dot_esubc;
            d_dot_esubc.x=d.x*e_sub_c.x;
            d_dot_esubc.y=d.y*e_sub_c.y;
            d_dot_esubc.z=d.z*e_sub_c.z;
            
            ddotesubc=d_dot_esubc.x+d_dot_esubc.y+d_dot_esubc.z;
            
            
            // ray d_dot_d;
            d_dot_d.x=d.x*d.x;
            d_dot_d.y=d.x*d.y;
            d_dot_d.z=d.x*d.z;
            
            ddotd=d_dot_d.x+d_dot_d.y+d_dot_d.z;
            
            //   ray e_sub_cdote_sub_c;
            e_sub_cdote_sub_c.x=e_sub_c.x*e_sub_c.x;
            e_sub_cdote_sub_c.y=e_sub_c.y*e_sub_c.y;
            e_sub_cdote_sub_c.z=e_sub_c.z*e_sub_c.z;
            
            e_cdote_c_minus_rsquare=e_sub_cdote_sub_c.x+e_sub_cdote_sub_c.y+e_sub_cdote_sub_c.z-rsquare;
            
            
            disc=(ddotesubc*ddotesubc)-(ddotd*e_cdote_c_minus_rsquare);
            if(disc<0)
            {
                
            }
            else
            {   t3=(-(ddotesubc)+sqrt(disc))/ddotd;
                t4=(-(ddotesubc)-sqrt(disc))/ddotd;
                if(t3>t4)
                {
                    if(t4<min)
                    {
                        min=t4;
                    }
                    pt2.x=e.x+(t4*d.x);
                    pt2.y=e.y+(t4*d.y);
                    pt2.z=e.z+(t4*d.z);
                    
                }
                else
                {
                    if(t3<min)
                    {
                        min=t3;
                    }
                    pt2.x=e.x+(t3*d.x);
                    pt2.y=e.y+(t3*d.y);
                    pt2.z=e.z+(t3*d.z);
                    
                }
                n.x=pt2.x-c.x;
                n.y=pt2.y-c.y;
                n.z=pt2.z-c.z;
                
                deno=sqrt((n.x*n.x)+(n.y*n.y)+(n.z*n.z));
                n.x=n.x/deno;
                n.y=n.y/deno;
                n.z=n.z/deno;
                
                ray n_dot_l;
                n_dot_l.x=n.x*l.x;
                n_dot_l.y=n.y*l.y;
                n_dot_l.z=n.z*l.z;
                
                ndotl=n_dot_l.x+n_dot_l.y+n_dot_l.z;
                
                if(ndotl>0.0)
                {
                    max=ndotl;
                }
                else
                    max=0.0;
                l_ds.r=kd*i_ds.x*max;
                l_ds.g=kd*i_ds.y*max;
                l_ds.b=kd*i_ds.z*max;
                
                h.x=d.x+l.x;
                h.y=d.y+l.y;
                h.z=d.z+l.z;
                
                deno=sqrt((h.x*h.x)+(h.y*h.y)+(h.z*h.z));
                h.x=h.x/deno;
                h.y=h.y/deno;
                h.z=h.z/deno;
                
                
                
                
                
                ray n_dot_h;
                n_dot_h.x=n.x*h.x;
                n_dot_h.y=n.y*h.y;
                n_dot_h.z=n.z*h.z;
                
                ndoth=n_dot_h.x+n_dot_h.y+n_dot_h.z;
                
                if(ndoth>0.0)
                {
                    max=ndoth;
                }
                else
                {
                    max=0.0;
                }
                
                l_sp.r=ks*i_sp.x*pow(max,n_spec_power);
                l_sp.g=ks*i_sp.y*pow(max,n_spec_power);
                l_sp.b=ks*i_sp.z*pow(max,n_spec_power);
                
                
                
                
                
                
                
                
                
                /*l_shading.r=l_sp.r;
                 l_shading.g=l_sp.g;
                 l_shading.b=l_sp.b;*/
                l_shading.r=l_ds.r+l_ambient.r+l_sp.r;
                l_shading.g=l_ds.g+l_ambient.g+l_sp.g;
                l_shading.b=l_ds.b+l_ambient.b+l_sp.b;
                
                
                
                
                
                
                /* l_shading.r=l_ds.r+l_ambient.r;
                 l_shading.g=l_ds.g+l_ambient.g;
                 l_shading.b=l_ds.b+l_ambient.b;*/
                
                
                /*pixels[i][j].r=1.0*l_ambient.r;
                 pixels[i][j].g=0.0*l_ambient.g;
                 pixels[i][j].b=0.0*l_ambient.b;
                 
                 
                 pixels[i][j].r=1.0;
                 pixels[i][j].g=0.0;
                 pixels[i][j].b=0.0; */
                
                pixels[i][j].r=0.0*l_shading.r;
                pixels[i][j].g=0.0*l_shading.g;
                pixels[i][j].b=1.0*l_shading.b;
            }
            
            
            A=a.x-b.x;
            B=a.y-b.y;
            C=a.z-b.z;
            D=a.x-ct1.x;
            E=a.y-ct1.y;
            F=a.z-ct1.z;
            G=d.x;
            H=d.y;
            I=d.z;
            J=a.x-e.x;
            
            e.y=i;
            e.z=j;
            K=a.y-e.y;
            L=a.z-e.z;
            M=(A*((E*I)-(H*F)))+(B*((G*F)-(D*I)))+(C*((D*H)-(E*G)));
            beta=((J*((E*I)-(H*F)))+(K*((G*F)-(D*I)))+(L*((D*H)-(E*G))))/M;
            gamma=((I*((A*K)-(J*B)))+(H*((J*C)-(A*L)))+(G*((B*L)-(K*C))))/M;
            t5=-((F*((A*K)-(J*B)))+(E*((J*C)-(A*L)))+(D*((B*L)-(K*C))))/M;
            if(t5<min)
            {
                min=t5;
            }
            if((beta>0)&&(gamma>0)&&((beta+gamma)<1))
            {
                
                A1.x=b.x-a.x;
                A1.y=b.y-a.y;
                A1.z=b.z-a.z;
                
                B1.x=b.x-ct1.x;
                B1.y=b.y-ct1.y;
                B1.z=b.z-ct1.z;
                
                N1.x=(A1.y*B1.z)-(A1.z*B1.y);
                N1.y=(A1.z*B1.x)-(A1.x*B1.z);
                N1.z=(A1.x*B1.y)-(A1.y*B1.x);
                
                
                deno=sqrt((N1.x*N1.x)+(N1.y*N1.y)+(N1.z*N1.z));
                N1.x=N1.x/deno;
                N1.y=N1.y/deno;
                N1.z=N1.z/deno;
                
                
                C1.x=a.x-b.x;
                C1.y=a.y-b.y;
                C1.z=a.z-b.z;
                
                D1.x=b.x-ct1.x;
                D1.y=b.y-ct1.y;
                D1.z=b.z-ct1.z;
                
                N2.x=(C1.y*D1.z)-(C1.z*D1.y);
                N2.y=(C1.z*D1.x)-(C1.x*D1.z);
                N2.z=(C1.x*D1.y)-(C1.y*D1.x);
                
                
                
                
                deno=sqrt((N2.x*N2.x)+(N2.y*N2.y)+(N2.z*N2.z));
                N2.x=N2.x/deno;
                N2.y=N2.y/deno;
                N2.z=N2.z/deno;
                
                E1.x=ct2.x-a.x;
                E1.y=ct2.y-a.y;
                E1.z=ct2.z-a.z;
                
                
                deno=sqrt((E1.x*E1.x)+(E1.y*E1.y)+(E1.z*E1.z));
                E1.x=E1.x/deno;
                E1.y=E1.y/deno;
                E1.z=E1.z/deno;
                
                costheta1=(E1.x*N1.x)+(E1.y*N1.y)+(E1.z*N1.z);
                costheta2=(E1.x*N2.x)+(E1.y*N2.y)+(E1.z*N2.z);
                
                if(costheta1>costheta2)
                {
                    N=N2;
                }
                else
                    N=N1;
                
                L_diffuse.r=l.x-(e.x+(t5*d.x));
                L_diffuse.g=l.y-(e.y+(t5*d.y));
                L_diffuse.b=l.z-(e.z+(t5*d.z));
                
                deno=sqrt((L_diffuse.r*L_diffuse.r)+(L_diffuse.g*L_diffuse.g)+(L_diffuse.b*L_diffuse.b));
                L_diffuse.r=L_diffuse.r/deno;
                L_diffuse.g=L_diffuse.g/deno;
                L_diffuse.b=L_diffuse.b/deno;
                
                
                NdotLdiffuse=N.x*L_diffuse.r+N.y*L_diffuse.g+N.z*L_diffuse.b;
                if(NdotLdiffuse>0)
                {
                    max=NdotLdiffuse;
                }
                else
                    max=0.0;
                
                
                l_ds.r=kd*i_ds.x*max;
                l_ds.g=kd*i_ds.y*max;
                l_ds.b=kd*i_ds.z*max;
                
                
                
                
                
                h.x=d.x+l.x;
                h.y=d.y+l.y;
                h.z=d.z+l.z;
                
                deno=sqrt((h.x*h.x)+(h.y*h.y)+(h.z*h.z));
                h.x=h.x/deno;
                h.y=h.y/deno;
                h.z=h.z/deno;
                
                
                
                
                
                ray n_dot_h;
                n_dot_h.x=N.x*h.x;
                n_dot_h.y=N.y*h.y;
                n_dot_h.z=N.z*h.z;
                
                ndoth=n_dot_h.x+n_dot_h.y+n_dot_h.z;
                
                if(ndoth>0.0)
                {
                    max=ndoth;
                }
                else
                {
                    max=0.0;
                }
                
                l_sp.r=ks*i_sp.x*pow(max,n_spec_power);
                l_sp.g=ks*i_sp.y*pow(max,n_spec_power);
                l_sp.b=ks*i_sp.z*pow(max,n_spec_power);
                
                
                
                
                
                
                
                
                
                
                l_shading.r=l_ds.r+l_ambient.r+l_sp.r;
                l_shading.g=l_ds.g+l_ambient.g+l_sp.g;
                l_shading.b=l_ds.b+l_ambient.b+l_sp.b;
                
                
                
                
                
                
                pixels[i][j].r=0.0*l_shading.r;
                pixels[i][j].g=1.0*l_shading.g;
                pixels[i][j].b=0.0*l_shading.b;
                
            }
            D=a.x-ct2.x;
            E=a.y-ct2.y;
            F=a.z-ct2.z;
            
            e.y=i;
            e.z=j;
            K=a.y-e.y;
            L=a.z-e.z;
            M=(A*((E*I)-(H*F)))+(B*((G*F)-(D*I)))+(C*((D*H)-(E*G)));
            beta=((J*((E*I)-(H*F)))+(K*((G*F)-(D*I)))+(L*((D*H)-(E*G))))/M;
            gamma=((I*((A*K)-(J*B)))+(H*((J*C)-(A*L)))+(G*((B*L)-(K*C))))/M;
            t6=-((F*((A*K)-(J*B)))+(E*((J*C)-(A*L)))+(D*((B*L)-(K*C))))/M;
            if(t6<min)
            {
                min=t6;
            }
            if((beta>0)&&(gamma>0)&&((beta+gamma)<1))
            {
                
                
                
                A1.x=b.x-a.x;
                A1.y=b.y-a.y;
                A1.z=b.z-a.z;
                
                B1.x=b.x-ct2.x;
                B1.y=b.y-ct2.y;
                B1.z=b.z-ct2.z;
                
                N1.x=(A1.y*B1.z)-(A1.z*B1.y);
                N1.y=(A1.z*B1.x)-(A1.x*B1.z);
                N1.z=(A1.x*B1.y)-(A1.y*B1.x);
                
                
                deno=sqrt((N1.x*N1.x)+(N1.y*N1.y)+(N1.z*N1.z));
                N1.x=N1.x/deno;
                N1.y=N1.y/deno;
                N1.z=N1.z/deno;
                
                
                C1.x=a.x-b.x;
                C1.y=a.y-b.y;
                C1.z=a.z-b.z;
                
                D1.x=b.x-ct2.x;
                D1.y=b.y-ct2.y;
                D1.z=b.z-ct2.z;
                
                N2.x=(C1.y*D1.z)-(C1.z*D1.y);
                N2.y=(C1.z*D1.x)-(C1.x*D1.z);
                N2.z=(C1.x*D1.y)-(C1.y*D1.x);
                
                
                
                
                deno=sqrt((N2.x*N2.x)+(N2.y*N2.y)+(N2.z*N2.z));
                N2.x=N2.x/deno;
                N2.y=N2.y/deno;
                N2.z=N2.z/deno;
                
                E1.x=ct1.x-ct2.x;
                E1.y=ct1.y-ct2.y;
                E1.z=ct1.z-ct2.z;
                
                
                deno=sqrt((E1.x*E1.x)+(E1.y*E1.y)+(E1.z*E1.z));
                E1.x=E1.x/deno;
                E1.y=E1.y/deno;
                E1.z=E1.z/deno;
                
                costheta1=E1.x*N1.x+E1.y*N1.y+E1.z*N1.z;
                costheta2=E1.x*N2.x+E1.y*N2.y+E1.z*N2.z;
                
                if(costheta1>costheta2)
                {
                    N=N2;
                }
                else
                    N=N1;
                
                L_diffuse.r=l.x-(e.x+(t6*d.x));
                L_diffuse.g=l.y-(e.y+(t6*d.y));
                L_diffuse.b=l.z-(e.z+(t6*d.z));
                
                deno=sqrt((L_diffuse.r*L_diffuse.r)+(L_diffuse.g*L_diffuse.g)+(L_diffuse.b*L_diffuse.b));
                L_diffuse.r=L_diffuse.r/deno;
                L_diffuse.g=L_diffuse.g/deno;
                L_diffuse.b=L_diffuse.b/deno;
                
                
                NdotLdiffuse=N.x*L_diffuse.r+N.y*L_diffuse.g+N.z*L_diffuse.b;
                if(NdotLdiffuse>0)
                {
                    max=NdotLdiffuse;
                }
                else
                    max=0.0;
                
                
                l_ds.r=kd*i_ds.x*max;
                l_ds.g=kd*i_ds.y*max;
                l_ds.b=kd*i_ds.z*max;
                
                
                
                h.x=d.x+l.x;
                h.y=d.y+l.y;
                h.z=d.z+l.z;
                
                deno=sqrt((h.x*h.x)+(h.y*h.y)+(h.z*h.z));
                h.x=h.x/deno;
                h.y=h.y/deno;
                h.z=h.z/deno;
                
                
                
                
                
                ray n_dot_h;
                n_dot_h.x=N.x*h.x;
                n_dot_h.y=N.y*h.y;
                n_dot_h.z=N.z*h.z;
                
                ndoth=n_dot_h.x+n_dot_h.y+n_dot_h.z;
                
                if(ndoth>0.0)
                {
                    max=ndoth;
                }
                else
                {
                    max=0.0;
                }
                
                l_sp.r=ks*i_sp.x*pow(max,n_spec_power);
                l_sp.g=ks*i_sp.y*pow(max,n_spec_power);
                l_sp.b=ks*i_sp.z*pow(max,n_spec_power);
                
                
                
                
                
                
                
                
                
                /*l_shading.r=l_sp.r;
                 l_shading.g=l_sp.g;
                 l_shading.b=l_sp.b;*/
                l_shading.r=l_ds.r+l_ambient.r+l_sp.r;
                l_shading.g=l_ds.g+l_ambient.g+l_sp.g;
                l_shading.b=l_ds.b+l_ambient.b+l_sp.b;
                
                
                
                /*
                 l_shading.r=l_ds.r;
                 l_shading.g=l_ds.g;
                 l_shading.b=l_ds.b;
                 
                 
                 l_shading.r=l_ds.r+l_ambient.r;
                 l_shading.g=l_ds.g+l_ambient.g;
                 l_shading.b=l_ds.b+l_ambient.b;*/
                
                
                pixels[i][j].r=0.0*l_shading.r;
                pixels[i][j].g=1.0*l_shading.g;
                pixels[i][j].b=0.0*l_shading.b;
                
                
                
                
            }
            A=a.x-ct2.x;
            B=a.y-ct2.y;
            C=a.z-ct2.z;
            D=a.x-ct1.x;
            E=a.y-ct1.y;
            F=a.z-ct1.z;
            G=d.x;
            H=d.y;
            I=d.z;
            J=a.x-e.x;
            
            e.y=i;
            e.z=j;
            K=a.y-e.y;
            L=a.z-e.z;
            M=(A*((E*I)-(H*F)))+(B*((G*F)-(D*I)))+(C*((D*H)-(E*G)));
            beta=((J*((E*I)-(H*F)))+(K*((G*F)-(D*I)))+(L*((D*H)-(E*G))))/M;
            gamma=((I*((A*K)-(J*B)))+(H*((J*C)-(A*L)))+(G*((B*L)-(K*C))))/M;
            t7=-((F*((A*K)-(J*B)))+(E*((J*C)-(A*L)))+(D*((B*L)-(K*C))))/M;
            if(t7<min)
            {
                min=t7;
            }
            if((beta>0)&&(gamma>0)&&((beta+gamma)<1))
            {
                A1.x=ct2.x-a.x;
                A1.y=ct2.y-a.y;
                A1.z=ct2.z-a.z;
                
                B1.x=ct2.x-ct1.x;
                B1.y=ct2.y-ct1.y;
                B1.z=ct2.z-ct1.z;
                
                N1.x=(A1.y*B1.z)-(A1.z*B1.y);
                N1.y=(A1.z*B1.x)-(A1.x*B1.z);
                N1.z=(A1.x*B1.y)-(A1.y*B1.x);
                
                
                deno=sqrt((N1.x*N1.x)+(N1.y*N1.y)+(N1.z*N1.z));
                N1.x=N1.x/deno;
                N1.y=N1.y/deno;
                N1.z=N1.z/deno;
                
                
                C1.x=a.x-ct2.x;
                C1.y=a.y-ct2.y;
                C1.z=a.z-ct2.z;
                
                D1.x=ct2.x-ct1.x;
                D1.y=ct2.y-ct1.y;
                D1.z=ct2.z-ct1.z;
                
                N2.x=(C1.y*D1.z)-(C1.z*D1.y);
                N2.y=(C1.z*D1.x)-(C1.x*D1.z);
                N2.z=(C1.x*D1.y)-(C1.y*D1.x);
                
                
                
                
                deno=sqrt((N2.x*N2.x)+(N2.y*N2.y)+(N2.z*N2.z));
                N2.x=N2.x/deno;
                N2.y=N2.y/deno;
                N2.z=N2.z/deno;
                
                E1.x=b.x-ct1.x;
                E1.y=b.y-ct1.y;
                E1.z=b.z-ct1.z;
                
                
                deno=sqrt((E1.x*E1.x)+(E1.y*E1.y)+(E1.z*E1.z));
                E1.x=E1.x/deno;
                E1.y=E1.y/deno;
                E1.z=E1.z/deno;
                
                costheta1=E1.x*N1.x+E1.y*N1.y+E1.z*N1.z;
                costheta2=E1.x*N2.x+E1.y*N2.y+E1.z*N2.z;
                
                if(costheta1>costheta2)
                {
                    N=N2;
                }
                else
                    N=N1;
                
                L_diffuse.r=l.x-(e.x+(t7*d.x));
                L_diffuse.g=l.y-(e.y+(t7*d.y));
                L_diffuse.b=l.z-(e.z+(t7*d.z));
                
                deno=sqrt((L_diffuse.r*L_diffuse.r)+(L_diffuse.g*L_diffuse.g)+(L_diffuse.b*L_diffuse.b));
                L_diffuse.r=L_diffuse.r/deno;
                L_diffuse.g=L_diffuse.g/deno;
                L_diffuse.b=L_diffuse.b/deno;
                
                
                NdotLdiffuse=N.x*L_diffuse.r+N.y*L_diffuse.g+N.z*L_diffuse.b;
                if(NdotLdiffuse>0)
                {
                    max=NdotLdiffuse;
                }
                else
                    max=0.0;
                
                
                l_ds.r=kd*i_ds.x*max;
                l_ds.g=kd*i_ds.y*max;
                l_ds.b=kd*i_ds.z*max;
                
                
                
                
                h.x=d.x+l.x;
                h.y=d.y+l.y;
                h.z=d.z+l.z;
                
                deno=sqrt((h.x*h.x)+(h.y*h.y)+(h.z*h.z));
                h.x=h.x/deno;
                h.y=h.y/deno;
                h.z=h.z/deno;
                
                
                
                
                
                ray n_dot_h;
                n_dot_h.x=N.x*h.x;
                n_dot_h.y=N.y*h.y;
                n_dot_h.z=N.z*h.z;
                
                ndoth=n_dot_h.x+n_dot_h.y+n_dot_h.z;
                
                if(ndoth>0.0)
                {
                    max=ndoth;
                }
                else
                {
                    max=0.0;
                }
                
                l_sp.r=ks*i_sp.x*pow(max,n_spec_power);
                l_sp.g=ks*i_sp.y*pow(max,n_spec_power);
                l_sp.b=ks*i_sp.z*pow(max,n_spec_power);
                
                
                
                
                
                
                
                
                
                /*l_shading.r=l_sp.r;
                 l_shading.g=l_sp.g;
                 l_shading.b=l_sp.b;*/
                l_shading.r=l_ds.r+l_ambient.r+l_sp.r;
                l_shading.g=l_ds.g+l_ambient.g+l_sp.g;
                l_shading.b=l_ds.b+l_ambient.b+l_sp.b;
                
                
                
                
                
                
                /* l_shading.r=l_ds.r+l_ambient.r;
                 l_shading.g=l_ds.g+l_ambient.g;
                 l_shading.b=l_ds.b+l_ambient.b; */
                
                
                pixels[i][j].r=0.0*l_shading.r;
                pixels[i][j].g=1.0*l_shading.g;
                pixels[i][j].b=0.0*l_shading.b;
                
            }
            A=ct1.x-ct2.x;
            B=ct1.y-ct2.y;
            C=ct1.z-ct2.z;
            D=ct1.x-b.x;
            E=ct1.y-b.y;
            F=ct1.z-b.z;
            G=d.x;
            H=d.y;
            I=d.z;
            J=ct1.x-e.x;
            
            e.y=i;
            e.z=j;
            K=ct1.y-e.y;
            L=ct1.z-e.z;
            M=(A*((E*I)-(H*F)))+(B*((G*F)-(D*I)))+(C*((D*H)-(E*G)));
            beta=((J*((E*I)-(H*F)))+(K*((G*F)-(D*I)))+(L*((D*H)-(E*G))))/M;
            gamma=((I*((A*K)-(J*B)))+(H*((J*C)-(A*L)))+(G*((B*L)-(K*C))))/M;
            t8=-((F*((A*K)-(J*B)))+(E*((J*C)-(A*L)))+(D*((B*L)-(K*C))))/M;
            if(t8<min)
            {
                min=t8;
            }
            if((beta>0)&&(gamma>0)&&((beta+gamma)<1))
            {
                A1.x=ct2.x-ct1.x;
                A1.y=ct2.y-ct1.y;
                A1.z=ct2.z-ct1.z;
                
                B1.x=ct2.x-b.x;
                B1.y=ct2.y-b.y;
                B1.z=ct2.z-b.z;
                
                N1.x=(A1.y*B1.z)-(A1.z*B1.y);
                N1.y=(A1.z*B1.x)-(A1.x*B1.z);
                N1.z=(A1.x*B1.y)-(A1.y*B1.x);
                
                
                deno=sqrt((N1.x*N1.x)+(N1.y*N1.y)+(N1.z*N1.z));
                N1.x=N1.x/deno;
                N1.y=N1.y/deno;
                N1.z=N1.z/deno;
                
                
                C1.x=ct1.x-ct2.x;
                C1.y=ct1.y-ct2.y;
                C1.z=ct1.z-ct2.z;
                
                D1.x=ct2.x-b.x;
                D1.y=ct2.y-b.y;
                D1.z=ct2.z-b.z;
                
                N2.x=(C1.y*D1.z)-(C1.z*D1.y);
                N2.y=(C1.z*D1.x)-(C1.x*D1.z);
                N2.z=(C1.x*D1.y)-(C1.y*D1.x);
                
                
                
                
                deno=sqrt((N2.x*N2.x)+(N2.y*N2.y)+(N2.z*N2.z));
                N2.x=N2.x/deno;
                N2.y=N2.y/deno;
                N2.z=N2.z/deno;
                
                E1.x=a.x-b.x;
                E1.y=a.y-b.y;
                E1.z=a.z-b.z;
                
                
                deno=sqrt((E1.x*E1.x)+(E1.y*E1.y)+(E1.z*E1.z));
                E1.x=E1.x/deno;
                E1.y=E1.y/deno;
                E1.z=E1.z/deno;
                
                costheta1=E1.x*N1.x+E1.y*N1.y+E1.z*N1.z;
                costheta2=E1.x*N2.x+E1.y*N2.y+E1.z*N2.z;
                
                if(costheta1>costheta2)
                {
                    N=N2;
                }
                else
                    N=N1;
                
                
                L_diffuse.r=l.x-(e.x+(t8*d.x));
                L_diffuse.g=l.y-(e.y+(t8*d.y));
                L_diffuse.b=l.z-(e.z+(t8*d.z));
                
                deno=sqrt((L_diffuse.r*L_diffuse.r)+(L_diffuse.g*L_diffuse.g)+(L_diffuse.b*L_diffuse.b));
                L_diffuse.r=L_diffuse.r/deno;
                L_diffuse.g=L_diffuse.g/deno;
                L_diffuse.b=L_diffuse.b/deno;
                
                
                NdotLdiffuse=N.x*L_diffuse.r+N.y*L_diffuse.g+N.z*L_diffuse.b;
                if(NdotLdiffuse>0)
                {
                    max=NdotLdiffuse;
                }
                else
                    max=0.0;
                
                
                l_ds.r=kd*i_ds.x*max;
                l_ds.g=kd*i_ds.y*max;
                l_ds.b=kd*i_ds.z*max;
                
                
                
                
                h.x=d.x+l.x;
                h.y=d.y+l.y;
                h.z=d.z+l.z;
                
                deno=sqrt((h.x*h.x)+(h.y*h.y)+(h.z*h.z));
                h.x=h.x/deno;
                h.y=h.y/deno;
                h.z=h.z/deno;
                
                
                
                
                
                ray n_dot_h;
                n_dot_h.x=N.x*h.x;
                n_dot_h.y=N.y*h.y;
                n_dot_h.z=N.z*h.z;
                
                ndoth=n_dot_h.x+n_dot_h.y+n_dot_h.z;
                
                if(ndoth>0.0)
                {
                    max=ndoth;
                }
                else
                {
                    max=0.0;
                }
                
                l_sp.r=ks*i_sp.x*pow(max,n_spec_power);
                l_sp.g=ks*i_sp.y*pow(max,n_spec_power);
                l_sp.b=ks*i_sp.z*pow(max,n_spec_power);
                
                
                
                
                
                
                
                
                
                /*l_shading.r=l_sp.r;
                 l_shading.g=l_sp.g;
                 l_shading.b=l_sp.b;*/
                l_shading.r=l_ds.r+l_ambient.r+l_sp.r;
                l_shading.g=l_ds.g+l_ambient.g+l_sp.g;
                l_shading.b=l_ds.b+l_ambient.b+l_sp.b;
                
                
                
                /*
                 l_shading.r=l_ds.r+l_ambient.r;
                 l_shading.g=l_ds.g+l_ambient.g;
                 l_shading.b=l_ds.b+l_ambient.b; */
                
                
                pixels[i][j].r=0.0*l_shading.r;
                pixels[i][j].g=1.0*l_shading.g;
                pixels[i][j].b=0.0*l_shading.b;
                
                if(min<40000)
                {
                    if((min==t1)||(min==t2))
                    {
                        pixels[i][j]=v1;
                    }
                    else if((min==t3)||(min==t4))
                    {
                        pixels[i][j]=v2;
                    }
                    else if((min==t5)||(min==t6)||(min==t7)||(min==t8))
                    {
                        pixels[i][j]=v3;
                    }
                    
                }
                
                
                
                
            }
            
            
            
            
        }
        
        
        
    }
    
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawPixels(512,512,GL_RGB,GL_FLOAT,pixels);
    
    
    glutSwapBuffers();
    //    pos+=delta;
    //    if(pos>=1.0 || pos<=-1.0)
    //    { delta=-delta; }
    //    glutPostRedisplay();
    // pos=+delta;
    // glutPostRedisplay();
}

int main(int argc,char**argv)
{
    glutInit( &argc, argv );
    glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE);
    glutInitWindowSize(512,512);
    glutCreateWindow("glut1");
    glutDisplayFunc( display );
    glutMainLoop();
    return 0;
}



