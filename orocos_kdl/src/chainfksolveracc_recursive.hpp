#include "chainfksolver.hpp"
#include <kdl/chain.hpp>
#include <kdl/chainfksolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainjnttojacdotsolver.hpp>

namespace KDL
{
    class ChainFkSolverAcc_recursive : public ChainFkSolverAcc
    {
    public:
        ChainFkSolverAcc_recursive(const KDL::Chain &chain)
            : m_chain(chain),
              jacdot_solver(chain),
              jac_solver(chain)
        {
        }

        int JntToCart(const JntArrayAcc &q_in, FrameAcc &out, int segmentNr = -1) override
        {
            const KDL::JntArray &kdl_q = q_in.q;
            const KDL::JntArray &kdl_qdot = q_in.qdot;
            const KDL::JntArray &kdl_qdotdot = q_in.qdotdot;

            KDL::JntArrayVel kdl_vel(kdl_q, kdl_qdot);

            int res;

            // Get jacdot*qdot
            KDL::Twist jdot_qdot;
            res = jacdot_solver.JntToJacDot(kdl_vel, jdot_qdot, segmentNr);
            if (0 != res)
            {
                return res;
            }

            // Get the Jacobian
            KDL::Jacobian jac(m_chain.getNrOfJoints());
            res = jac_solver.JntToJac(kdl_q, jac, segmentNr);
            if (0 != res)
            {
                return res;
            }

            // EEF acceleration is jac*qdotdot + jacdot*qdot
            KDL::Twist kdl_eef_accel;
            KDL::MultiplyJacobian(jac, kdl_qdotdot, kdl_eef_accel);
            kdl_eef_accel = kdl_eef_accel + jdot_qdot;

            out.p.dv = kdl_eef_accel.vel;
            out.M.dw = kdl_eef_accel.rot;

            return 0;
        }

        int JntToCart(const JntArrayAcc& q_in, std::vector<FrameAcc>& out,int segmentNr=-1) override
        {
            for (auto i = 0; i < m_chain.getNrOfSegments(); ++i)
            {
                if (int res = JntToCart(q_in, out[i], i))
                {
                    return res;
                }
            }

            return 0;
        }
    
        void updateInternalDataStructures() override
        {

        }

        const KDL::Chain &m_chain;
        KDL::ChainJntToJacDotSolver jacdot_solver;
        KDL::ChainJntToJacSolver jac_solver;
    };
} // namespace KDL
