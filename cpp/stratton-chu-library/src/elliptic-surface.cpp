#include "stratton-chu/elliptic-surface.hpp"
#include <iostream>

// Constructor of half of Elliptic surface r = r0 + α * x + β * y + [α×β] / |[α×β]| * c*(1 - x²/a² - y²/b²)^0.5
EllipticSurface::EllipticSurface(const Position& r0, const Vector& alpha, const Vector& beta,
                                   const double a, const double b, const double c) :
    m_r0(r0), m_alpha(alpha), m_beta(beta), m_a(a), m_b(b), m_c(c)
{
    m_n = m_alpha % m_beta;
    m_n /= m_n.norm();
}

bool EllipticSurface::BelongsToSurface (const Vector2D& pos) {
    if (sqr(pos[0])/sqr(m_a)+sqr(pos[1])/sqr(m_b)>1.0) {return false;}
    else {return true;}
}

Position EllipticSurface::point(const Vector2D& pos) const
{
    Position result;

    if (sqr(pos[0])/sqr(m_a)+sqr(pos[1])/sqr(m_b)<=1.0) {
        result = m_n;
        result *= m_c*sqrt(1.0 - sqr(pos[0])/sqr(m_a) - sqr(pos[1])/sqr(m_b));
        result += m_r0 + m_alpha * pos[0] + m_beta * pos[1];
        // std::cout << result[0] << " " << result[1] << " " << result[2] << std::endl;
    }
    else {
        // std::cout << "BAD POINT APPEARED" << std::endl;
        double newr=m_b*sqrt(1.0 - sqr(pos[0])/sqr(m_a));
        // std::cout << newr << " " << pos[1]  <<" "<< pos[0] <<std::endl;
        // double newz=sqrt(newr*newr-sqr(pos[1]));
        
        // double newpos1=pos[1]/fabs(pos[1])*sqrt(m_b*m_b*(1.0 - sqr(pos[0])/sqr(m_a))-newz*newz);
        // double newz=sqr(newr)/sqrt(sqr(newr)+sqr(pos[1]));
        // double newpos1=newz*pos[1]/newr;
        double newpos1=pos[1]/fabs(pos[1])*sqrt(9*sqr(newr)*sqr(pos[1])/(sqr(newr)+9*sqr(pos[1])));

        // std::cout << pos[0] << " " << pos[1] << " " << newz << std::endl;
        result = m_n;
        // std::cout << result[0] << " " << result[1] << " " << result[2] << std::endl;
        result *= m_c*sqrt(1.0 - sqr(pos[0])/sqr(m_a) - sqr(newpos1)/sqr(m_b));
        // std::cout << sqrt(1.0 - sqr(pos[0])/sqr(m_a) - sqr(newpos1)/sqr(m_b)) <<std::endl;
        // result *= newz;
        // std::cout << result[0] << " " << result[1] << " " << result[2] << std::endl;
        // result += m_r0 + m_alpha * pos[0] + m_beta * pos[1];
        result += m_r0 + m_alpha * pos[0] + m_beta * newpos1;
        // std::cout << result[0] << " " << result[1] << " " << result[2] << std::endl;

        if (sqr(pos[0])/sqr(m_a)+sqr(newpos1)/sqr(m_b)==1.0) {
            std::cout << "lol" << std::endl;
        }
        
        // result[0]=0.0;
        // result[1]=0.0;
        // result[2]=0.0;
    }   
  
    return result;
}

Vector2D EllipticSurface::point2d(const Vector& pos) const
{
    Vector2D result;
    Vector support=pos-m_r0;
    result[0] = support*m_alpha;
    result[1] = support*m_beta;
    if (abs(result[0])>m_a || abs(result[1])>m_b) {std::cout << "the point is outside of ellipsoid" << std::endl; }
    return result;
}


Vector EllipticSurface::tau1(const Vector2D& pos) const
{
    Vector r_x;
    if (sqr(pos[0])/sqr(m_a)+sqr(pos[1])/sqr(m_b)>1.0) {
        // std::cout << "BAD POINT INPUT tau1" << std::endl;
        // double newpos1=pos[1]/fabs(pos[1])*m_b*sqrt(1.0 - sqr(pos[0])/sqr(m_a));

        // std::cout << pos[0] << " " << pos[1] << " " << newz << std::endl;
        // result = m_n;
        // std::cout << result[0] << " " << result[1] << " " << result[2] << std::endl;
        // // result *= m_c*sqrt(1.0 - sqr(pos[0])/sqr(m_a) - sqr(newpos1)/sqr(m_b));
        // result *= m_c*newz;
        // std::cout << result[0] << " " << result[1] << " " << result[2] << std::endl;
        // result += m_r0 + m_alpha * pos[0] + m_beta * pos[1];
        // result = m_r0 + m_alpha * pos[0] + m_beta * newpos1;
        double newr=m_b*sqrt(1.0 - sqr(pos[0])/sqr(m_a));
        double newpos1=pos[1]/fabs(pos[1])*sqrt(9*sqr(newr)*sqr(pos[1])/(sqr(newr)+9*sqr(pos[1])));
        r_x = m_n;
        r_x *= -m_c*pos[0] / sqr(m_a) /sqrt((1.0 - sqr(pos[0])/sqr(m_a) - sqr(newpos1)/sqr(m_b)));
        r_x += m_alpha;
        // r_x[0]=0.0;
        // r_x[1]=0.0;
        // r_x[2]=0.0;
    }
    else if (sqr(pos[0])/sqr(m_a)+sqr(pos[1])/sqr(m_b)==1.0) {
        std::cout << "edge point (review the normal vector)" << std::endl;
        r_x = m_alpha*(-pos[1]/sqr(m_b))+m_beta*(pos[0]/sqr(m_a));
    } else {
        r_x = m_n;
        r_x *= -m_c*pos[0] / sqr(m_a) /sqrt((1.0 - sqr(pos[0])/sqr(m_a) - sqr(pos[1])/sqr(m_b)));
        r_x += m_alpha;
    }    
    
    return r_x;
}

Vector EllipticSurface::tau2(const Vector2D& pos) const
{
    Vector r_y;

    if (sqr(pos[0])/sqr(m_a)+sqr(pos[1])/sqr(m_b)>1.0) {
        // std::cout << "BAD POINT INPUT" << std::endl;
        double newr=m_b*sqrt(1.0 - sqr(pos[0])/sqr(m_a));
        double newpos1=pos[1]/fabs(pos[1])*sqrt(9*sqr(newr)*sqr(pos[1])/(sqr(newr)+9*sqr(pos[1])));
        r_y = m_n;
        r_y *= -m_c*newpos1 / sqr(m_b) /sqrt((1.0 - sqr(pos[0])/sqr(m_a) - sqr(newpos1)/sqr(m_b)));
        r_y += m_beta;
        // r_y[0]=0.0;
        // r_y[1]=0.0;
        // r_y[2]=0.0;
    }
    else if (sqr(pos[0])/sqr(m_a)+sqr(pos[1])/sqr(m_b)==1.0) {
        std::cout << "edge point (review the normal vector)" << std::endl;
        r_y = m_n;
    } else {
        r_y = m_n;
        r_y *= -m_c*pos[1] / sqr(m_b) /sqrt((1.0 - sqr(pos[0])/sqr(m_a) - sqr(pos[1])/sqr(m_b)));
        r_y += m_beta;
    }  

    return r_y;
}

Vector EllipticSurface::dS_over_dxdy(const Vector2D& pos) const {
    Vector norm;
    // norm=(tau1(pos) % tau2(pos));
    // std::cout << "norm bad " << norm[0] << " " << norm[1] << " " << norm[2] <<  std::endl;
    norm=-(tau1(pos) % tau2(pos));
    // std::cout << "pos " << pos[0] <<" "<< pos[1]<<" norm " << norm[0] << " " << norm[1] << " " << norm[2] <<  std::endl;
    return norm;
}