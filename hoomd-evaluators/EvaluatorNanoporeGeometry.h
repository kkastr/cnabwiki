// Copyright (c) 2009-2017 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: jglaser

#ifndef __EVALUATOR_NANOPORE_GEOMETRY_H__
#define __EVALUATOR_NANOPORE_GEOMETRY_H__

#ifndef NVCC
#include <string>
#endif

#include <math.h>
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

/*! \file EvaluatorExternalElectricField.h
    \brief Defines the external potential evaluator to induce a periodic ordered phase
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define DEVICE __device__
#else
#define DEVICE
#endif

//! Class for constructing the geometry of a nanopore out of wca potentials.

class EvaluatorNanoporeGeometry
    {
    public:

        //! type of parameters this external potential accepts
        typedef Scalar4 param_type;
        typedef Scalar3 field_type;

        //! Constructs the constraint evaluator
        /*! \param X position of particle
            \param box box dimensions
            \param params per-type parameters of external potential
        */
        DEVICE EvaluatorNanoporeGeometry(Scalar3 X, const BoxDim& box, const param_type& _params, const field_type& field)
            : m_pos(X),
              m_box(box),
              m_field(field), lj1(_params.x), lj2(_params.y), rpore(_params.z), tpore(_params.w)
            {
            }

        //! External Periodic doesn't need diameters
        DEVICE static bool needsDiameter() { return false; }
        //! Accept the optional diameter value
        /*! \param di Diameter of particle i
        */
        DEVICE void setDiameter(Scalar di) { }

        //! External Periodic doesn't need charges
        DEVICE static bool needsCharge() { return false; }
        //! Accept the optional diameter value
        /*! \param qi Charge of particle i
        */
        DEVICE void setCharge(Scalar qi) { }

        //! Declares additional virial cotribututions are needed for the external field
        /*! No contribution
        */
        DEVICE static bool requestFieldVirialTerm() { return true; }

        //! Evaluate the force, energy and virial
        /*! \param F force vector
            \param energy value of the energy
            \param virial array of six scalars for the upper triangular virial tensor
        */
        DEVICE void evalForceEnergyAndVirial(Scalar3& F, Scalar& energy, Scalar* virial)
            {


            Scalar cx = m_pos.x - m_field.x;
            Scalar cy = m_pos.y - m_field.y;
            Scalar cz = m_pos.z - m_field.z;
            // Scalar r = sqrt(cx*cx + cy*cy +cz*cz);
            Scalar rxy = sqrt(cx*cx + cy*cy);



            // WCA parameters


            Scalar rtop ;
            Scalar rbot ;

            Scalar rtop2;
            Scalar rbot2;

            Scalar rtop2inv;
            Scalar rbot2inv;

            Scalar WCAforcemag = Scalar(0.0);

            // WCA Cutoff
            Scalar rcut = 1.12246204831;

            // Inter-pore parameters
            Scalar theta;
            Scalar htpore = tpore*0.5;
            Scalar zmax = htpore + rcut;
            Scalar rxy2 = (rxy-rpore)*(rxy-rpore);
            Scalar rxy2inv = Scalar(1.0)/rxy2;
            Scalar rxy6inv = rxy2inv * rxy2inv * rxy2inv;
            Scalar phi;
            Scalar R2inv;
            Scalar R6inv;
            Scalar Fmag;
            Scalar R2;




            // Force initialization
            F.x = 0;
            F.y = 0;
            F.z = 0;
            Scalar Ftheta;

            if (rxy <= rpore-rcut){
                F.x = 0;
                F.y = 0;
                F.z = 0;



            } else if (abs(cz) >= zmax ){
                    F.x = 0;
                    F.y = 0;
                    F.z = 0;
                    } else if (rxy >= rpore && abs(cz) >= htpore) {

                    if (cz > 0) {
                        rtop = cz - (htpore);
                        rtop2 =rtop*rtop;
                        rtop2inv = Scalar(1.0)/rtop2;
                        Scalar rtop6inv= rtop2inv*rtop2inv*rtop2inv;
                        WCAforcemag = rtop*rtop2inv * rtop6inv * (Scalar(12.0)*lj1*rtop6inv - Scalar(6.0)*lj2);
                        //pair_eng = (r6inv * (lj1*r6inv - lj2) + epsilon);
                        F.x = F.y = 0;
                        F.z = WCAforcemag;
                        // F.z = 0;
                    // }
                    } else {
                        rbot = cz + (htpore);
                        rbot2 =rbot*rbot;
                        rbot2inv = Scalar(1.0)/rbot2;
                        Scalar rbot6inv= rbot2inv*rbot2inv*rbot2inv;
                        WCAforcemag = rbot*rbot2inv * rbot6inv * (Scalar(12.0)*lj1*rbot6inv - Scalar(6.0)*lj2);
                        //pair_eng = (r6inv * (lj1*r6inv - lj2) + epsilon);
                        F.x = F.y = 0;
                        F.z =  WCAforcemag;
                        // F.z = 0;
                    }

                } else if (abs(cz) <= htpore  && rxy > rpore - rcut){
                        theta = atan2(cy,cx);

                         Ftheta = (rxy-rpore)*rxy2inv * rxy6inv * (Scalar(12.0)*lj1*rxy6inv - Scalar(6.0)*lj2);
                         F.z = 0 ;
                         F.x = Ftheta*cos(theta);
                         F.y = Ftheta*sin(theta);



                } else if (abs(cz) >= htpore && rxy > rpore -rcut) {

                theta = atan2(cy,cx);
                Scalar rxycorner = sqrt((cx-rpore*cos(theta))*(cx-rpore*cos(theta)) + (cy-rpore*sin(theta))*(cy-rpore*sin(theta)));
                R2 = (rxycorner)*(rxycorner) + ((abs(cz)-htpore))*((abs(cz)-htpore));
                //Scalar theta2 = acos(((cz-tpore))/sqrt(R2));

                    if (R2 <= rcut*rcut){
                        R2inv = Scalar(1.0)/R2;
                        R6inv = R2inv*R2inv*R2inv;
                        Fmag  = sqrt(R2)*R2inv * R6inv * (Scalar(12.0)*lj1*R6inv - Scalar(6.0)*lj2);



                        if (cz >= 0) {
                            phi = atan2(cz - htpore,rxycorner);
                            F.x = -Fmag*cos(phi)*cos(theta);
                            F.y = -Fmag*sin(theta)*cos(phi);
                            F.z = Fmag*sin(phi);
                        } else {
                            phi = atan2(cz + htpore,rxycorner);
                            F.x = -Fmag*cos(phi)*cos(theta);
                            F.y = -Fmag*sin(theta)*cos(phi);
                            F.z = Fmag*sin(phi);

                        }

                    } else { F.x = F.y= F.z= 0 ;}



                }


            

            virial[0] = F.x*m_pos.x;
            virial[1] = F.x*m_pos.y;
            virial[2] = F.x*m_pos.z;
            virial[3] = F.y*m_pos.y;
            virial[4] = F.y*m_pos.z;
            virial[5] = F.z*m_pos.z;   
            }


        #ifndef NVCC
        //! Get the name of this potential
        /*! \returns The potential name. Must be short and all lowercase, as this is the name energies will be logged as
            via analyze.log.
        */
        static std::string getName()
            {
            return std::string("nanopore");
            }
        #endif

    protected:
        Scalar3 m_pos;                //!< particle position
        BoxDim m_box;                 //!< box dimensions
        Scalar m_qi;                  //!< particle charge
        Scalar3 m_field;              //!< the field vector
        Scalar lj1;                   //!< lj1 parameter extracted from the params passed to the constructor
        Scalar lj2;                   //!< lj2 parameter extracted from the params passed to the constructor
        Scalar rpore;                 //!< Radius of Pore
        Scalar tpore;                 //!< Thickness of Pore
   };


#endif // __EVALUATOR_EXTERNAL_LAMELLAR_H__
