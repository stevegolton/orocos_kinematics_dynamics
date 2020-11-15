#include <kdl/chain.hpp>
#include <cmath>
#include <chainfksolveracc_recursive.hpp>
#include <iostream>

int main(void)
{
    using namespace KDL;

    const double PI_2 = M_PI_2;

    Chain puma560;
    puma560.addSegment(Segment());
    puma560.addSegment(Segment(Joint(Joint::RotZ),
                               Frame::DH(0.0, PI_2, 0.0, 0.0),
                               RigidBodyInertia(0, Vector::Zero(), RotationalInertia(0, 0.35, 0, 0, 0, 0))));
    puma560.addSegment(Segment(Joint(Joint::RotZ),
                               Frame::DH(0.4318, 0.0, 0.0, 0.0),
                               RigidBodyInertia(17.4, Vector(-.3638, .006, .2275), RotationalInertia(0.13, 0.524, 0.539, 0, 0, 0))));
    puma560.addSegment(Segment());
    puma560.addSegment(Segment(Joint(Joint::RotZ),
                               Frame::DH(0.0203, -PI_2, 0.15005, 0.0),
                               RigidBodyInertia(4.8, Vector(-.0203, -.0141, .070), RotationalInertia(0.066, 0.086, 0.0125, 0, 0, 0))));
    puma560.addSegment(Segment(Joint(Joint::RotZ),
                               Frame::DH(0.0, PI_2, 0.4318, 0.0),
                               RigidBodyInertia(0.82, Vector(0, .019, 0), RotationalInertia(1.8e-3, 1.3e-3, 1.8e-3, 0, 0, 0))));
    puma560.addSegment(Segment());
    puma560.addSegment(Segment());
    puma560.addSegment(Segment(Joint(Joint::RotZ),
                               Frame::DH(0.0, -PI_2, 0.0, 0.0),
                               RigidBodyInertia(0.34, Vector::Zero(), RotationalInertia(.3e-3, .4e-3, .3e-3, 0, 0, 0))));
    puma560.addSegment(Segment(Joint(Joint::RotZ),
                               Frame::DH(0.0, 0.0, 0.0, 0.0),
                               RigidBodyInertia(0.09, Vector(0, 0, .032), RotationalInertia(.15e-3, 0.15e-3, .04e-3, 0, 0, 0))));
    puma560.addSegment(Segment());


    ChainFkSolverAcc_recursive solver(puma560);

    JntArrayAcc q(puma560.getNrOfJoints());
    q.qdotdot(0) = 1.25;
    q.qdotdot(1) = 1.25;
    q.qdotdot(2) = 1.25;
    q.qdotdot(3) = 1.25;
    q.qdotdot(4) = 1.25;
    q.qdotdot(5) = 1.25;
    
    std::vector<FrameAcc> f(puma560.getNrOfSegments());

    int res = solver.JntToCart(q, f);

    if (res)
    {
        std::cerr << "err" << '\n';
    }

    std::cout << "Acc = " << f[puma560.getNrOfSegments() - 1].p.dv.x() << '\n';

    return 0;
}