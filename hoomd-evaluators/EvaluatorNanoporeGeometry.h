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
        DEVICE static bool requestFieldVirialTerm() { return false; }

        //! Evaluate the force, energy and virial
        /*! \param F force vector
            \param energy value of the energy
            \param virial array of six scalars for the upper triangular virial tensor
        */
        DEVICE void evalForceEnergyAndVirial(Scalar3& F, Scalar& energy, Scalar* virial)
            {
            
            // WCA Cutoff 2^(1/6)
            
            Scalar r_cut = 1.12246204831;
            
            Scalar cx = m_pos.x - m_field.x;
            Scalar cy = m_pos.y - m_field.y;
            Scalar cz = m_pos.z - m_field.z;
            
            Scalar half_tpore = tpore*0.5;
            Scalar zmax = half_tpore + r_cut;


            Scalar sgn_cz = cz/abs(cz);

            Scalar dz= sgn_cz * (abs(cz)-half_tpore);
            Scalar dr;
            

            // WCA force init
            Scalar f_wca_interior = Scalar(0.0);
            Scalar f_wca_exterior = Scalar(0.0);
            Scalar f_wca_corner = Scalar(0.0);
            Scalar dz7inv;
            Scalar dz13inv;

            Scalar dr7inv;
            Scalar dr13inv;

            Scalar drho7inv;
            Scalar drho13inv;

           
            
            
            // inside the pore 
            Scalar r_xy = fast::sqrt(cx*cx + cy*cy);
            Scalar r_interior = rpore - r_cut;
           

            // corner of the pore
            Scalar x_corner;
            Scalar y_corner;
            Scalar r_corner;
            Scalar rho_corner;

            Scalar phi;
            Scalar theta_corner;

            Scalar theta = atan2(cy,cx);
            
            x_corner = cx - rpore*fast::cos(theta);
            y_corner = cy - rpore*fast::sin(theta);
           
            r_corner = fast::sqrt(x_corner*x_corner + y_corner*y_corner);
            
            rho_corner = fast::sqrt(x_corner*x_corner + y_corner*y_corner + dz*dz);

            // Initialize forces
            F.x = F.y = F.z =0;

            if (r_xy<r_interior || abs(cz)>=zmax ){
                F.x = F.y = F.z =0;
            } else if (r_xy >= rpore && abs(cz)<zmax && abs(cz) >= half_tpore){

                dz7inv = Scalar(1.0)/pow(dz,7);
                dz13inv = Scalar(1.0)/pow(dz,13);
                
                f_wca_exterior = Scalar(12.0)*lj1*dz13inv - Scalar(6.0)*lj2*dz7inv;

                F.x = F.y =0;
                F.z = f_wca_exterior;
            } else if (r_xy >= r_interior && abs(cz)<= half_tpore) {

                
                dr = r_xy - rpore;

                dr7inv = Scalar(1.0)/pow(dr,7);
                dr13inv = Scalar(1.0)/pow(dr,13);

                f_wca_interior = Scalar(12.0)*lj1*dr13inv - Scalar(6.0)*lj2*dr7inv;


                F.x = f_wca_interior*fast::cos(theta);
                F.y = f_wca_interior*fast::sin(theta);
                F.z = 0;
            } else if (r_xy>r_interior && abs(cz)>half_tpore && rho_corner<=r_cut){
                

                theta_corner = atan2(y_corner,x_corner);
                phi = atan2(dz,r_corner);
                drho7inv = Scalar(1.0)/pow(rho_corner,7);
                drho13inv = Scalar(1.0)/pow(rho_corner,13);

                f_wca_corner = Scalar(12.0)*lj1*drho13inv - Scalar(6.0)*lj2*drho7inv;
                
                
                F.x = f_wca_corner * fast::cos(theta_corner)*fast::cos(phi);
                F.y = f_wca_corner * fast::sin(theta_corner)*fast::cos(phi);
                F.z = f_wca_corner * fast::sin(phi);
            } 



        





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
