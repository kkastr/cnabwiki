// Copyright (c) 2009-2017 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.



#ifndef __EVALUATOR_NANOPORE_FRICTION_H__
#define __EVALUATOR_NANOPORE_FRICTION_H__

#ifndef NVCC
#include <string>
#endif

#include <math.h>
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"


// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define DEVICE __device__
#else
#define DEVICE
#endif

//! Class for simulating friction on large particles from the pore

class EvaluatorNanoporeFriction
    {
    public:

        //! type of parameters this external potential accepts
       //!typedef struct param{} param_type;
        typedef Scalar2 param_type;
        typedef Scalar3 field_type;

        //! Constructs the constraint evaluator
        /*! \param X position of particle
            \param box box dimensions
            \param params per-type parameters of external potential
        */
        DEVICE EvaluatorNanoporeFriction(Scalar3 X, const BoxDim& box, const param_type& params, const field_type& field)
            : m_pos(X),
              m_box(box), cf(field), f_params(params)
            {
            }

        
        DEVICE static bool needsDiameter() { return true; }
        //! Accept the optional diameter value
        /*! \param di Diameter of particle i
        */
        DEVICE void setDiameter(Scalar di) { m_di = di; }

       
        DEVICE static bool needsCharge() { return false; }

        /*! \param qi Charge of particle i
        */
        DEVICE void setCharge(Scalar qi) {}

        //! Declares additional virial cotribututions are needed for the external field
        /*! No contribution
        */
        DEVICE static bool requestFieldVirialTerm() { return false; }

        //! Evaluate the force, energy and virial
        /*! \param F force vector
            \param energy value of the energy
            \param virial array of six scalars for the upper triangular virial tensor
        */
        DEVICE void evalForceEnergyAndVirial(Scalar3& F, Scalar& energy, Scalar* virial)
            {

            // cx,cy,cz are the coordinates of the center of the nanopore


            Scalar rx = m_pos.x - cf.x;
            Scalar ry = m_pos.y - cf.y;
            Scalar rz = m_pos.z - cf.z;

            

            Scalar epsilon = f_params.x;
            Scalar alpha =  f_params.y;

            
            if (m_di > Scalar(1.0)){       

                Scalar dr = fast::sqrt(rx*rx +ry*ry +rz*rz);
                Scalar U = epsilon*((dr*dr)/(dr*dr +Scalar(1.0)/alpha) -1);  
                Scalar Fi = -epsilon * ((Scalar(2.0) * alpha*dr)/pow(alpha*dr*dr + 1,2));     
                F.x = Fi * rx/dr;
                F.y = Fi * ry/dr;
                F.z = Fi * rz/dr;
                energy = U;
            }
            else {

                F.x = 0;
                F.y = 0;
                F.z = 0;
                energy=0;
            }

      
            }

        #ifndef NVCC
        //! Get the name of this potential
        /*! \returns The potential name. Short and all lowercase.
        */
        static std::string getName()
            {
            return std::string("pore_friction");
            }
        #endif

    protected:
        Scalar3 m_pos;                //!< particle position
        BoxDim m_box;                 //!< box dimensions
        Scalar m_di;                  //!< particle diameter
        Scalar3 cf;                  //!< contains friction center
        Scalar2 f_params;             //!< contains alpha,epsilon,tpore
        

   };


#endif 
