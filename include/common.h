//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_COMMON_H
#define XROTOR_NOGRAPHICS_CPP_COMMON_H

#include <fstream>
#include <cstdio>
#include <cstdarg>
#include <string>
using namespace std;

namespace common {
    const int IX = 100;                     /**< max number of radial prop stations */
    const int IPX = IX + 1;
    const int NPARX = 12;                   /**< max number of stored cases */
    const int ICASX = 100;                  /**< number of case parameters stored */
    const int IWX = 200;                    /**< dimension of work arrays */

    const int NAX = 20;                     /**< max number of aerodynamic sections defined */
    const int NDX = 14;                     /**< number of aerodynamic parameters defined for each section */

    const int IQ = IX + 5;
    const int JX = (IQ * 3) / 2 + 1;

    const double PI = 3.14159265358979323846264338327950288;

    struct context {
        double  Q           [IQ][IQ];       /**< working array */

        double  MCRIT;                      /**< critical Mach number */

        bool    CONV,                       /**< true if Converged solution exists */
                GREEK,                      /**< true if Unrecognized command */
                TERSE,                      /**< true if Terse output (no radial distributions) */
                VRTX,                       /**<  */
                FAST,                       /**< true if Graded Momentum, Potential otherwise */
                FREE,                       /**< true if Free wake */
                DUCT,                       /**< true if duct is present */
                LSTRUCT,                    /**< */
                DEST,                       /**< */
                DESP,                       /**< */
                LDESINI,                    /**< */
                LOPRINI,                    /**< */
                LROTOR,                     /**< */
                LVNORM,                     /**< */
                LPWRVAR,                    /**< */
                STALL       [IX],           /**< true for each station which is stalled */
                LEGEND;                     /**< */

        string  SAVIL,                      /**< Disk output save filename */
                FNAME,                      /**< Generic filename */
                NAME;                       /**< Case name */

        double  RHO,                        /**< Fluid density            (dimensioned) */
                RMU,                        /**< Fluid dynamic viscosity  (dimensioned) */
                VSO,                        /**< Fluid speed of sound     (dimensioned) */
                VEL,                        /**< Flight speed             (dimensioned) */
                RAD,                        /**< Rotor tip radius         (dimensioned) */
                GEE,                        /**< Earth's acceleration     (dimensioned) */
                ALT;                        /**< Altitude for fluid properties (km),  999.0 if not defined */

        int     II,                         /**< Number of radial stations on blade */
                IINF,                       /**< Number of radial stations on blade + outer domain */
                INCR,                       /**< Radial station increment for terminal output */
                NN,                         /**< Number of perturbation potential harmonic terms */
                NBLDS,                      /**< Number of blades */
                IXSPAC,                     /**< 1 = cosine r/R array stretching
                                             *   2 = sine stretching (less spacing near root)*/
                NITERD,                     /**< */
                NITERA;                     /**< */

        float   VERSION;                    /**< */
        double  DT;                         /**< */

        int     IAERO       [IX];           /**< */

        double  CH          [IX],           /**< Chord array */
                BETA        [IX],           /**< Twist angle array */
                BETA0       [IX],           /**< Static twist angle array (with zero structural loads) */
                T           [IX],           /**< Dummy radial coordinate array */
                DBETA,                      /**< Accumulated change in twist angle */
                XI          [IX],           /**< Radial coordinate array (r/R) */
                DXI         [IX],           /**< Radial coordinate increment at each station */
                XI0,                        /**< Blade root radial coordinate value */
                XITIP,                      /**< Blade tip radial coordinate value  (always = 1) */
                XINF,                       /**< Outer radial coordinate value where farfield BC is applied */
                XPITCH,                     /**< x/c location of pitch axis for loads calculations and plots */
                XV          [IX],           /**< */
                RAKE;                       /**< */

        int     NADD;                       /**< */
        double  RADD        [IX],           /**< */
                UADD        [IX],           /**< */
                VADD        [IX],           /**< */
                UADDR       [IX],           /**< */
                VADDR       [IX],           /**< */
                UBODY       [IX],           /**< Nacelle perturbation axial  velocity */
                VBODY       [IX],           /**< Nacelle perturbation radial velocity */
                URDUCT;                     /**< */

        double  CL          [IX],           /**< Local lift coefficient array */
                CD          [IX],           /**< Local drag coefficient array */
                CM          [IX],           /**< Local blade airfoil Cm */
                RE          [IX],           /**< Local Reynolds number array */
                EFFP        [IX],           /**< Local profile efficiency array */
                GAM         [IX],           /**< Local circulation array */
                DTII        [IX],           /**< */
                DPII        [IX],           /**< */
                DTVI        [IX],           /**< */
                DPVI        [IX],           /**< */
                DTWI        [IX],           /**< */
                DPWI        [IX];           /**< */

        // aero data quantities for each defined radial aerodynamic section
        int     NAERO;                      /**< Number of aerodynamic datasets defined (NAERO>=1) */
        double  XIAERO      [NAX],          /**< Radial station r/R where aero dataset is defined */
                AERODATA    [NDX][NAX];     /**< Aerodynamic definition of the blade section at XIAERO
                                             *   AERODATA[ 0][:] = A0 (angle of zero lift)
                                             *   AERODATA[ 2][:] = CLMAX (Max CL)
                                             *   AERODATA[ 3][:] = CLMIN (Min CL)
                                             *   AERODATA[ 4][:] = DCLDA (Incompressible 2-D lift curve slope)
                                             *   AERODATA[ 5][:] = DCLDA_STALL (2-D lift curve slope at stall)
                                             *   AERODATA[ 6][:] = DCL_STALL (CL increment, onset to full stall)
                                             *   AERODATA[ 7][:] = CDMIN (Minimum drag coefficient value)
                                             *   AERODATA[ 8][:] = CLDMIN (Lift at minimum drag value)
                                             *   AERODATA[ 9][:] = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
                                             *   AERODATA[10][:] = CMCON (Incompressible 2-D pitching moment)
                                             *   AERODATA[11][:] = REREF (reference Reynold's number)
                                             *   AERODATA[12][:] = REXP (Reynold's number exponent Cd~Re^REXP)
                                             *   AERODATA[13][:] = MCRIT (critical Mach #)
                                             *   AERODATA[14][:] = TOC (thickness/chord)
                                             */

        int     NCASE,                      /**< current number of saved operating cases */
                KCASE,                      /**< indicator for independent parameter of case sweep:
                                             *       0 = none
                                             *       1 = advance ratio
                                             *       2 = rpm
                                             *       3 = velocity
                                             *       4 = blade angle
                                             */
                IWTYP;                      /**< Type of induced velocity model emplyed currently
                                             *   1 = Graded Momentum,  2 = Potential Formulation */

        double  ADV,                        /**< Advance ratio */
                ADW,                        /**< Wake advance ratio */
                ADWFCTR,                    /**< */
                RMS,                        /**< Residual */
                RLX,                        /**< Under-relaxation factor */
                EFFINV,                     /**< 1 / Inviscid efficiency */
                TSPEC,                      /**< Specified thrust */
                PSPEC,                      /**< Specified power */
                QSPEC,                      /**< Specified torque */
                TTOT,                       /**< Rotor inviscid + viscous + nacelle thrust */
                PTOT,                       /**< Rotor inviscid + viscous + nacelle power */
                QTOT,                       /**< Rotor inviscid + viscous + nacelle torque  = PTOT*ADV */
                TINV,                       /**< Inviscid thrust */
                PINV,                       /**< Inviscid power */
                TWAK,                       /**< Inviscid + nacelle thrust */
                PWAK,                       /**< Inviscid + nacelle power */
                TVIS,                       /**< Viscous thrust */
                PVIS,                       /**< Viscous power */
                                            // TTOT  =  TVIS +  TWAK
                                            //       =  TVIS + (TINV + Tnacelle)
                GRESMX,                     /**< */
                FRESMX,                     /**< */
                ARESMX;                     /**< */

        double  TI_ADV,                     /**< Sensitivity of TINV to advance ratio */
                PI_ADV,                     /**< Sensitivity of PINV to advance ratio */
                TI_ADW,                     /**< Sensitivity of TINV to wake advance ratio */
                PI_ADW,                     /**< Sensitivity of PINV to wake advance ratio */
                TW_ADV,                     /**< Sensitivity of TWAK to advance ratio */
                PW_ADV,                     /**< Sensitivity of PWAK to advance ratio */
                TW_ADW,                     /**< Sensitivity of TWAK to wake advance ratio */
                PW_ADW,                     /**< Sensitivity of PWAK to wake advance ratio */
                TV_ADV,                     /**< Sensitivity of TVIS to advance ratio */
                PV_ADV,                     /**< Sensitivity of PVIS to advance ratio */
                TV_ADW,                     /**< Sensitivity of TVIS to wake advance ratio */
                PV_ADW,                     /**< Sensitivity of PVIS to wake advance ratio */
                TV_DBE,                     /**< Sensitivity of TVIS to blade angle change */
                PV_DBE,                     /**< Sensitivity of PVIS to blade angle change */
                TI_GAM      [IX],           /**< Sensitivity of TINV to radial circulation distribution */
                PI_GAM      [IX],           /**< Sensitivity of PINV to radial circulation distribution */
                TW_GAM      [IX],           /**< Sensitivity of TWAK to radial circulation distribution */
                PW_GAM      [IX],           /**< Sensitivity of PWAK to radial circulation distribution */
                TV_GAM      [IX],           /**< Sensitivity of TVIS to radial circulation distribution */
                PV_GAM      [IX];           /**< Sensitivity of PVIS to radial circulation distribution */

        double  W0          [IWX],          /**< Temporary working array */
                W1          [IWX],          /**< Temporary working array */
                W2          [IWX],          /**< Temporary working array */
                W3          [IWX],          /**< Temporary working array */
                W4          [IWX],          /**< Temporary working array */
                W5          [IWX],          /**< Temporary working array */
                W6          [IWX],          /**< Temporary working array */
                W7          [IWX],          /**< Temporary working array */
                W8          [IWX],          /**< Temporary working array */
                W9          [IWX];          /**< Temporary working array */

        int     NPWRVAR;                    /**< */

        double  VWAK        [IX],           /**< Equivalent-rotor tangential induced velocity array */
                VW_GAM      [IX][IX],       /**< VWAK-blade bound circulation  sensitivity matrix */
                VW_ADW      [IX],           /**< VWAK-wake advance ratio sensitivity matrix */
                VW_ADV      [IX],           /**< VWAK-advance ratio sensitivity matrix */
                VIND        [3][IX],        /**< Tangential induced velocity array */
                VIND_GAM    [3][IX][IX],    /**< VIND-blade bound circulation  sensitivity matrix*/
                VIND_ADW    [3][IX];        /**< VIND-wake advance ratio sensitivity matrix*/

        double  XW          [IX],           /**< Equivalent-rotor radial coordinate array */
                XW_GAM      [IX][IX],       /**< XW-blade bound circulation  sensitivity matrix */
                XW_ADW      [IX][IX],       /**< XW-wake advance ratio sensitivity matrix */
                XW_ADV      [IX][IX],       /**< XW-advance ratio sensitivity matrix */
                DWX         [IX],           /**< Equivalent-rotor radial coordinate increment array */
                DXW_GAM     [IX][IX],       /**< DXW-blade bound circulation  sensitivity matrix */
                DXW_ADW     [IX][IX],       /**< DXW-wake advance ratio sensitivity matrix */
                DXW_ADV     [IX][IX];       /**< DXW-advance ratio sensitivity matrix */

        double  DGAM        [IX],           /**< Newton update delta array for bound circulation */
                RES         [IQ],           /**< Newton residual array for bound circulation */
                DADV,                       /**< Newton update delta for advance ratio */
                DADW,                       /**< Newton update delta for wake advance ratio */
                DBET,                       /**< Newton update delta for blade angle */
                DEFF,                       /**< Newton update delta for 1 / inviscid efficiency */
                DQ          [IQ],           /**< Generic solution vector */
                DGAMOLD     [IX];           /**< Previous value of DGAM */

    };

    /**
     * Write a formatted string to an output stream.
     *
     * @param os        output stream
     * @param format    format string
     * @param ...       arguments to the format string
     */
    void fprintf(ostream &os, const char* format, ...) {
        va_list argv;
        va_start(argv, format);
        char buff[1024];
        vsprintf(buff, format, argv);
        va_end(argv);
        os << buff;
    }
}






#endif //XROTOR_NOGRAPHICS_CPP_COMMON_H
