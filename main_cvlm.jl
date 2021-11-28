#Libraries used
using YAML
using Parameters
using Interpolations
using LinearAlgebra
using ProgressBars
using Plots
using NPZ
using DifferentialEquations

"""
    read_inp(debug=0)

Read the input from the input file, creates all the global numbering variables used in later stages

Inputs:
	inp_file- Name of the input YAML file with .yaml extension
    debug   - Integer 0,1 to switch on or off debugging

Outputs:
	sref      - reference surface area
	bref      - reference span
	cref      - reference chord length
	aseq      - sequence of alpha to investigate of form [starting,ending,number of steps] of variable types [float,float,integer]
	vinf      - value of vinf
	ρ         - density
	alpha     - angle of attack
	nwing     - number of wings
	ntot_patch- total number of patches
	npatch_perwing
	          - patches per wing (used for global patches numbering)
	ntot_lat  - total number of lattices (used for global lattice numbering)
	nspan_perpatch
	          - number of spanwise strips per patch (used to loop over spanwise strips)
	nlat_perpatch
			  - number of total lattices per patch
	xyz_qc_patch
			  - quarter chord locations for each patch of form [(x,y,z),(start,end),global patch number]
	chord_patch
		      - chord lengths for the pacthes of form [(start,end),global patch number]
	twist_patch
			  - twists for the pacthces of the form [(start,end),gloabl patch number]
	 α_zl_patch
	 		  - zero lift alpha values of the form [(start,end),gloabl patch number]

"""
function read_inp(inp_file,debug=0)
        inp_data = YAML.load_file(inp_file)
         @info "Reading Input File"
        name = inp_data["name"]
        sref = inp_data["sref"]
        bref = inp_data["bref"]
        cref = inp_data["cref"]
        nwing= inp_data["nwing"]

        #create arrays required to store
        npatch_perwing = zeros(Int64,nwing)
         for i in 1:nwing
                str = "wing"*string(i)
                npatch_perwing[i] = inp_data[str]["npatch"]
        end

        ntot_patch = sum(npatch_perwing)
        #get some arrays to store root and tip info
        xyz_qc_patch = zeros(3,2,ntot_patch)
        chord_patch  = zeros(2,ntot_patch)
        twist_patch  = zeros(2,ntot_patch)
        α_zl_patch   = zeros(2,ntot_patch)

        nlat_perpatch = zeros(Int64,ntot_patch)
		nspan_perpatch = zeros(Int64,ntot_patch)
        ntot_lat   = 0
        counter    = 0
        for i in 1:nwing
               str = "wing"*string(i)
                for j in 1:npatch_perwing[i]
                        counter = counter + 1
						nspan_perpatch[counter]   = inp_data[str]["patch"*string(j)]["nspan"]
                        nlat_perpatch[counter]    = inp_data[str]["patch"*string(j)]["nlat"]*nspan_perpatch[counter]

                        # at root
                        xyz_qc_patch[:,1,counter] = inp_data[str]["patch"*string(j)]["root"]["xyz_qc"]
                        chord_patch[1,counter]    = inp_data[str]["patch"*string(j)]["root"]["chord"]
                        twist_patch[1,counter]    = inp_data[str]["patch"*string(j)]["root"]["twist"]*pi/180
                        α_zl_patch[1,counter]     = inp_data[str]["patch"*string(j)]["root"]["alpha_zl"]*pi/180

                        # at tip
                        xyz_qc_patch[:,2,counter] = inp_data[str]["patch"*string(j)]["tip"]["xyz_qc"]
                        chord_patch[2,counter]    = inp_data[str]["patch"*string(j)]["tip"]["chord"]
                        twist_patch[2,counter]    = inp_data[str]["patch"*string(j)]["tip"]["twist"]*pi/180
                        α_zl_patch[2,counter]     = inp_data[str]["patch"*string(j)]["tip"]["alpha_zl"]*pi/180

                        ntot_lat = ntot_lat + nlat_perpatch[counter]
                end
        end

		aseq = inp_data["aseq"]
		aseq = [aseq[1]*pi/180,aseq[2]*pi/180,aseq[3]]
		vinf = inp_data["vinf"]
		ρ    = inp_data["rho"]
		alpha= inp_data["alpha"]*pi/180


        if (debug==1)
                #checks

                #print some stuff if needed
                display("ntot_lat:"*string(ntot_lat))
                display(xyz_qc_patch)

        end

        return (sref,bref,cref,aseq,vinf,ρ,alpha,nwing,ntot_patch,npatch_perwing,ntot_lat,nspan_perpatch,nlat_perpatch,xyz_qc_patch,chord_patch,twist_patch,α_zl_patch)
end

"""
    geom_calc(nwing,npatch,npatch_perwing,ntot_lat,nlat_perpatch,xyz_qc_patch,chord_patch,twist_patch,α_zl_patch,Λ,debug)

Calculate the geometric propoerties for all the lattice elements

	Inputs:
		nwing      - number of wings
		ntot_patch - total number of patches
		npatch_perwing
				   - number of patches per wing, form : [Integer value]
		ntot_lat   - total numer of lattices
		nspan_perpatch
				   - number of spanwise strips per patch
		nlat_perpatch
				   - number of lattices per patch (uses global patch numbering)
		xyz_qc_patch
				   - quarter chord locations for each patch of form [(x,y,z),(start,end),global patch number]
		 chord_patch
				   - chord lengths for the pacthes of form [(start,end),global patch number]
		 twist_patch
				   - twists for the pacthces of the form [(start,end),gloabl patch number]
		α_zl_patch
	 			   - zero lift alpha values of the form [(start,end),gloabl patch numb
		debug      - debug flag to enable some debugging options

	Output:
		nspn       - number of span wise locations to plt local distributions like Cl
						: If more than one wing, each wings span wise locations are stored in that order
						  Example - Wing + Tail
						  nspn = span locations on wing + span locations on tail
		spn_map    - mapping for all the lattices to span locations
						: form - Array{Integers},size "ntot_lat" pointing to position in "spn_loc" array
		spn_loc    - locations along span (y-coordinates)
						: form - Array{Float,1}, size "nspn"
		sbar       - locations of starting point of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		ebar       - locations of ending point of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		mbar       - locations of middle point of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		nbar       - perpendicular unit vectors of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		cbar       - locations of control points or quarter chord points of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		dbar       - locations of middle points of horse-shoe lattice for all lattices at Treffts Plane
						: form - Array{Float,2}, size {3,ntot_lat}
		tbar       - unit vector for pointing the tail of horse shoe vortices
						: form - Array{Float,2}, size {3,ntot_lat}
		chord      - chord length values at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		twist      - twist values at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		α_zl       - alpha zero lift values at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		θ          - dihedral angle at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		Λ          - sweep angle at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		ds         - length of filament for all lattices
						: form - Array{Float,1}, size {ntot_lat}

"""
function geom_calc(nwing,ntot_patch,npatch_perwing,ntot_lat,nspan_perpatch,nlat_perpatch,xyz_qc_patch,chord_patch,twist_patch,α_zl_patch,debug=0)

        #create arrays for the lattices
        sbar = zeros(3,ntot_lat)        # starting of bv
        ebar = zeros(3,ntot_lat)        # ending of bv
        mbar = zeros(3,ntot_lat)        # mid of bv
        cbar = zeros(3,ntot_lat)        # ctrl poi of bv
        dbar = zeros(3,ntot_lat)        # trefftz plane mid poi
        nbar = zeros(3,ntot_lat)		# normal for patch/wing/lattice
        tbar = zeros(3,ntot_lat)        # tail direction
		univ = zeros(3,ntot_lat)        #used for calc nbar

        chord = zeros(ntot_lat)         # chord length
        twist = zeros(ntot_lat)			# twist
        α_zl  = zeros(ntot_lat)			# α_zero-lift
        θ     = zeros(ntot_lat)			# dihedral angle
		Λ     = zeros(ntot_lat)         # sweep angle
		ds    = zeros(ntot_lat)         # lattice length


		spn_map = zeros(Int64,ntot_lat) #Spanwise mapping


		nlatspn_perwing = zeros(Int64,nwing) #Span per wing

		# Loop over to find total span wise locations (Could be eliminated/improved)
		nspn  = 0
		ipatch=0
		for i in 1:nwing
			for j in 1:npatch_perwing[i]
				ipatch = ipatch + 1
				nlat_perspan = Int64(nlat_perpatch[ipatch]/nspan_perpatch[ipatch])
				for k in 1:nlat_perspan
					nspn = nspn + 1
				end
			end
		end


		spn_loc = zeros(nspn) #Spanwise locations

		#Initialise counters
        ilat  = 0
        ipatch= 0
		nspn  = 0


        for i in  tqdm(1:nwing,leave=false)
                for j in 1:npatch_perwing[i]
                        ipatch = ipatch + 1 #for global patch number

						# Interpolate chord, twist and α_zl values
                        itp_chord = LinearInterpolation(xyz_qc_patch[2,:,ipatch],chord_patch[:,ipatch])
                        itp_twist = LinearInterpolation(xyz_qc_patch[2,:,ipatch],twist_patch[:,ipatch])
                        itp_α_zl  = LinearInterpolation(xyz_qc_patch[2,:,ipatch],α_zl_patch[:,ipatch])

						# Get the starting and ending locations along the chord for each patch
						nlat_perspan = Int64(nlat_perpatch[ipatch]/nspan_perpatch[ipatch])
						a = copy(xyz_qc_patch[:,1,ipatch])
                        b = copy(xyz_qc_patch[:,2,ipatch])

						#
						# Computing starting root and tip locations for each span
						#
						a_span = zeros(3,nspan_perpatch[ipatch])
						b_span = zeros(3,nspan_perpatch[ipatch])

						#pvec is perpendicular to (b-a)
						pvec   = (b[1:2] - a[1:2])/norm(b[1:2]-a[1:2])
						pvec   = [pvec[2],-pvec[1]]

						#loop over strips along the chord wise
						ns     = nspan_perpatch[ipatch]
						for ispan in 1:ns
							a_span[:,ispan] = a + 0.25*chord_patch[1,ipatch]*(4*ispan-3-ns)/ns * [pvec[1],pvec[2],0.0]
							b_span[:,ispan] = b + 0.25*chord_patch[2,ipatch]*(4*ispan-3-ns)/ns * [pvec[1],pvec[2],0.0]
						end

						#Sanity check to see the pvec vector is in the right direction
						@info pvec
						if (dot(pvec,[1.0,0.0])<=0.0)
							println("Error in geom_calc span wise variables")
							println("Perp Vector Calc wrong")
							return
						end

						#Loop over the chord wise strips and lattices in each strip
						for ispan in 1:nspan_perpatch[ipatch]
						for k in 1:nlat_perspan

                                ilat = ilat + 1

                                a = copy(a_span[:,ispan])
                                b = copy(b_span[:,ispan])

								#Find start,end,mid,normal vec, dihedral, sweep, α_zl, ...
								sbar[:,ilat] =  a + (b - a)/nlat_perspan .*(k-1)
                                ebar[:,ilat] =  a + (b - a)/nlat_perspan .*(k)
                                mbar[:,ilat] = 0.5*(sbar[:,ilat]+ebar[:,ilat])
                                univ[:,ilat] = cross((ebar[:,ilat] - sbar[:,ilat]),[1.0,0.0,0.0])/norm(cross((ebar[:,ilat] - sbar[:,ilat]),[-1.0,0.0,0.0]))

                                θ[ilat]      = acos(nbar[3,ilat])
								Λ[ilat]      = acos(dot(b[1:2]-a[1:2],[0.0,1.0])/norm(b[1:2]-a[1:2]))
                                chord[ilat]  = itp_chord(mbar[2,ilat])
                                twist[ilat]  = itp_twist(mbar[2,ilat])
                                α_zl[ilat]   = itp_α_zl(mbar[2,ilat])
								ds[ilat]	 = norm(ebar[2:3,ilat]-sbar[2:3,ilat])#*cos(Λ[ilat])

                                cbar[:,ilat] = mbar[:,ilat] - [chord[ilat]/(2*nspan_perpatch[ipatch]),0.0,0.0]
                                dbar[:,ilat] = [-10^8,mbar[2,ilat],mbar[3,ilat]]
								tbar[:,ilat] = [1.0,0.0,0.0]   #Trailing edge vector
								rot          = [cos(twist[ilat]) 0.0 sin(twist[ilat]);0.0 1.0 0.0;-sin(twist[ilat]) 0.0 cos(twist[ilat])]
								nbar[:,ilat] = rot*univ[:,ilat]

								#To get the spn_loc locations only when using the first strip
								if (ispan==1)
									nspn = nspn + 1
									spn_loc[nspn] = mbar[2,ilat]
								end

								#Compute ispn_bef to locate where in spn_map this patch is starting
								ispn_bef = 0
								for ip in 1:ipatch-1
									ispn_bef = ispn_bef + Int64(nlat_perpatch[ip]/nspan_perpatch[ip])
								end
								spn_map[ilat] = ispn_bef + k

						end
						end
					#finding numer of lattices in spn wise per each wing
					nlatspn_perwing[i] = nlatspn_perwing[i] + nlat_perspan
				end
        end

		#Debug stuff
        if debug==1
                display("sbar")
                display(sbar)
                display("xyz_qc_patch")
                display(xyz_qc_patch)
        end

        return nspn,spn_map,spn_loc,sbar,ebar,mbar,nbar,cbar,dbar,tbar,chord,twist,α_zl,θ,Λ,ds
end
"""
    calc_vind_sinf(rbar,sbar,tbar)

    Function to influence for semi-infinte vortex element given starting,ending and tail direction

    Input :
		rbar - vector for pointing the location of evaluation
				: form - Array{Float,1}, size {3}
		sbar - vector for pointing the start of filament
				: form - Array{Float,1}, size {3}
		tbar - vector for pointing the trailing of vortex
				: form - Array{Float,1}, size {3}

    Output:
		q    - vector of influence at the given location
				: form - Array{Float,1}, size {3}

"""
function calc_q_sinf(rbar,sbar,tbar)

	(x1,y1,z1) = sbar
	(x,y,z)    = rbar

	avec = [x1-x,y1-y,z1-z]


	amag = norm(avec)


	cprd = cross(avec,tbar)
	dprd = dot(avec,tbar)

	num = cprd* (1.0-(dprd/(amag)))
	den = dot(cprd,cprd)

	#for the point coinciding with the element
	if (den<1e-6)
		q = [0.0,0.0,0.0]
	else
		q =  num/den
	end
	return q
end
"""
    calc_vind_finite(rbar,sbar,ebar)

    Function to find the three velocity componenets for downwash given the line segment and the point

	Input :
		rbar - vector for pointing the location of evaluation
				: form - Array{Float,1}, size {3}
		sbar - vector for pointing the start of filament
				: form - Array{Float,1}, size {3}
		ebar - vector for pointing the end of filament
				: form - Array{Float,1}, size {3}

    Output:
		q    - vector of influence at the given location
				: form - Array{Float,1}, size {3}
"""
function calc_q_finite(rbar,sbar,ebar)

	(x1,y1,z1) = sbar
	(x2,y2,z2) = ebar
	(x,y,z)    = rbar

	avec = [x1-x,y1-y,z1-z]
	bvec = [x2-x,y2-y,z2-z]

	amag = norm(avec)
	bmag = norm(bvec)

	cprd = cross(avec,bvec)
	dprd = dot(avec,bvec)

	num = cprd* (amag+bmag) * (1.0-(dprd/(amag*bmag)))
	den = dot(cprd,cprd)

	#for the point coinciding with the element
	if (den<1e-6)
		q = [0.0,0.0,0.0]
	else
		q =  num/den
	end
	return q
end
"""
    calc_AICs()

    Calculate the AIC matrices required
"""
function calc_AICs(ntot_lat,sbar,ebar,cbar,mbar,nbar,tbar,debug=0)

	pi_inv = 1.0/pi

        AIC  = zeros(ntot_lat,ntot_lat)

		#Other matrices for CVLM
		Gxgamma = zeros(ntot_lat,ntot_lat)
		Gygamma = zeros(ntot_lat,ntot_lat)
		Gzgamma = zeros(ntot_lat,ntot_lat)
		Hxgamma = zeros(ntot_lat,ntot_lat)
		Hygamma = zeros(ntot_lat,ntot_lat)
		Hzgamma = zeros(ntot_lat,ntot_lat)

        AICₘ = zeros(ntot_lat,ntot_lat)
        AICₜ = zeros(ntot_lat,ntot_lat)

	xyz_ref = zeros(3,1)
	ubar = ebar - sbar
	r_bk = cbar - (ones(ntot_lat,1)*xyz_ref')'
	r_bj = mbar - (ones(ntot_lat,1)*xyz_ref')'
	for i in 1:ntot_lat
		for j in 1:ntot_lat

			#Calculating AIC
			#bounded vortex contribution
			qbv     = calc_q_finite(cbar[:,i],sbar[:,j],ebar[:,j])
			qbv_comp= dot(qbv,nbar[:,i])
			#left leg contribution
			qlv     = -calc_q_sinf(cbar[:,i],sbar[:,j],tbar[:,j])
			qlv_comp= dot(qlv,nbar[:,i])
			#right leg contribution
			qrv     = calc_q_sinf(cbar[:,i],ebar[:,j],tbar[:,j])
			qrv_comp= dot(qrv,nbar[:,i])

			AIC[i,j] = 0.25*pi_inv*(qbv_comp+qlv_comp+qrv_comp)
			qbv .= 0.0
			qbv_comp = 0.0
			qlv .= 0.0
			qlv_comp = 0.0
			qrv .= 0.0
			qrv_comp = 0.0

			#Calculating AICₘ
			#bounded vortex contribution
			qbv     = calc_q_finite(mbar[:,i],sbar[:,j],ebar[:,j])
			qbv_comp= dot(qbv,-nbar[:,i])
			#left leg contribution
			qlv     = -calc_q_sinf(mbar[:,i],sbar[:,j],tbar[:,j])
			qlv_comp= dot(qlv,-nbar[:,i])
			#right leg contribution
			qrv     = calc_q_sinf(mbar[:,i],ebar[:,j],tbar[:,j])
			qrv_comp= dot(qrv,-nbar[:,i])

			h         = 0.25*pi_inv*(qbv+qlv+qrv)
			AICₘ[i,j] = 0.25*pi_inv*(qbv_comp+qlv_comp+qrv_comp)

			Gxgamma[i,j] = dot(cross(h,ubar[:,j]),[1.0 0.0 0.0])
			Gygamma[i,j] = dot(cross(h,ubar[:,j]),[0.0 1.0 0.0])
			Gzgamma[i,j] = dot(cross(h,ubar[:,j]),[0.0 0.0 1.0])
			Hxgamma[i,j] = dot(cross(r_bj[:,j],cross(h,ubar[:,j])),[1.0 0.0 0.0])
			Hygamma[i,j] = dot(cross(r_bj[:,j],cross(h,ubar[:,j])),[0.0 1.0 0.0])
			Hzgamma[i,j] = dot(cross(r_bj[:,j],cross(h,ubar[:,j])),[0.0 0.0 1.0])


			qbv .= 0.0
			qbv_comp = 0.0
			qlv .= 0.0
			qlv_comp = 0.0
			qrv .= 0.0
			qrv_comp = 0.0

			#Calculating AICₜ
			#coordinate transformation for calculation purposes
			c0 = cbar[:,i]
			s0 = sbar[:,j]
			e0 = ebar[:,j]

			c0[1] = s0[1] = e0[1] = 0.0

			#bounded vortex contribution
			qbv     = 0.0
			#left leg contribution
			qlv     = -calc_q_sinf(c0,s0,tbar[:,j])
			qlv_comp= dot(qlv,nbar[:,i])
			#right leg contribution
			qrv     = calc_q_sinf(c0,e0,tbar[:,j])
			qrv_comp= dot(qrv,nbar[:,i])

			AICₜ[i,j] = 0.25*pi_inv*(qbv_comp+2.0*qlv_comp+2.0*qrv_comp)
		end
	end

	Gxv = zeros(ntot_lat,3)
	Gyv = zeros(ntot_lat,3)
	Gzv = zeros(ntot_lat,3)
	Hxv = zeros(ntot_lat,3)
	Hyv = zeros(ntot_lat,3)
	Hzv = zeros(ntot_lat,3)

	Gxomega = zeros(ntot_lat,3)
	Gyomega = zeros(ntot_lat,3)
	Gzomega = zeros(ntot_lat,3)
	Hxomega = zeros(ntot_lat,3)
	Hyomega = zeros(ntot_lat,3)
	Hzomega = zeros(ntot_lat,3)

	crsbar  = zeros(3,ntot_lat)

	for i in 1:ntot_lat
		crsbar[:,i] = cross(r_bk[:,i],nbar[:,i])

		Gxv[i,:] = cross(ubar[:,i],[1.0;0.0;0.0])
		Gyv[i,:] = cross(ubar[:,i],[0.0;1.0;0.0])
		Gzv[i,:] = cross(ubar[:,i],[0.0;0.0;1.0])

		Gxomega[i,:] = cross(r_bj[:,i],Gxv[i,:])
		Gyomega[i,:] = cross(r_bj[:,i],Gyv[i,:])
		Gzomega[i,:] = cross(r_bj[:,i],Gzv[i,:])

		Hxv[i,:] = cross(-ubar[:,i],cross(r_bj[:,i],[1.0;0.0;0.0]))
		Hyv[i,:] = cross(-ubar[:,i],cross(r_bj[:,i],[0.0;1.0;0.0]))
		Hzv[i,:] = cross(-ubar[:,i],cross(r_bj[:,i],[0.0;0.0;1.0]))

		Hxomega[i,:] = cross(-r_bj[:,i],cross(ubar[:,i],cross(r_bj[:,i],[1.0;0.0;0.0])))
		Hyomega[i,:] = cross(-r_bj[:,i],cross(ubar[:,i],cross(r_bj[:,i],[0.0;1.0;0.0])))
		Hzomega[i,:] = cross(-r_bj[:,i],cross(ubar[:,i],cross(r_bj[:,i],[0.0;0.0;1.0])))
	end

	J = AIC\[nbar' crsbar']

	P_x = J'*(-[Gxv Gxomega]+Gxgamma*J)
	P_y = J'*(-[Gyv Gyomega]+Gygamma*J)
	P_z = J'*(-[Gzv Gzomega]+Gzgamma*J)

	Q_x = J'*(-[Hxv Hxomega]+Hxgamma*J)
	Q_y = J'*(-[Hyv Hyomega]+Hygamma*J)
	Q_z = J'*(-[Hzv Hzomega]+Hzgamma*J)

    return AIC,AICₘ,AICₜ,J,P_x,P_y,P_z,Q_x,Q_y,Q_z,Gxgamma,Gygamma,
	Gzgamma,Hxgamma,Hygamma,Hzgamma,Gxv,Gyv,Gzv,Hxv,Hyv,Hzv,Gxomega,
	Gyomega,Gzomega,Hxomega,Hyomega,Hzomega,crsbar
end
"""
    calc_motion()

	Calculating the 6DOF of motion

	Input:

	Output:
"""
function calc_motion(Px,Py,Pz,Qx,Qy,Qz)

	# time span
	t_begin = 0
	t_end   = 1.0
	tsteps  = 11
	tspan = (t_begin,t_end)

	##########################################
	# initial conditions
	##########################################
	u0 = 100
	v0 = 0
	w0 = 3

	p0 = 0
	q0 = 0
	r0 = 0

	phi0 = 0
	theta0 = deg2rad(2)
	psi0 = 0

	x0 = 0
	y0 = 0
	z0 = -500


	###############
	# Let's get it running with CVLM in the loop - then, we will need to come up with a trim routine here for steady level flight. I've done this before, so I can add it.
	##############

	# inital condition vector [tranlational velocities | angular velocities | euler angles | cartesian position ]
	xst0 = [u0;v0;w0;p0;q0;r0;phi0;theta0;psi0;x0;y0;z0]

	#############################################
	# parameter definitions
	#############################################
	m = 10
	g = 9.81
	rho = 1.225
	T = 0.2 # thrust (will need to determine from induced drag at trim?)
	# Inertia constants for calculating GAMMA, C1-C9
	Ix = 1
	Iy = 1
	Iz = 1
	Ixz = 0.2
	# constants to deal with inertial coupling
	GAMMA = Ix*Iz-Ixz^2
	c1 = ((Iy-Iz)*Iz-Ixz^2)/GAMMA
	c2 = ((Ix-Iy+Iz)*Ixz)/GAMMA
	c3 =  Iz/GAMMA
	c4 = Ixz/GAMMA
	c5 = (Iz-Ix)/Iy
	c6 = Ixz/Iy
	c7 = 1/Iy
	c8 = ((Ix-Iy)*Ix-Ixz^2)/GAMMA
	c9 = Ix/GAMMA

	# rotating mass along longitudinal axis (spinning propellor or turbomachinery) - these will probably always be zero for using
	Ip = 0  # rotating inertia
	OmegaP = 0  # rpm of rotating machinery



	##################################################
	# Define parameter vector (would be more convienient as a tuple)
	##################################################
	pq = [m,g,rho,c1,c2,c3,c4,c5,c6,c7,c8,c9,T,Ip,OmegaP]


	## ODE Function
	function calc_dx!(dxst,xst,pq,t)

	    #X = [u,v,w,p,q,r,phi,theta,psi,x,y,z]

	    # parameters:
	    m,g,rho,c1,c2,c3,c4,c5,c6,c7,c8,c9,T,Ip,OmegaP = pq
	    # current state values:
	    u,v,w,p,q,r,phi,theta,psi,x,y,z = xst

	    # Aerodynamic Force/Moment Calculation (without controls {for now})
	    X = [u;v;w;p;q;r]

	    Fx = rho*X'*Px*X
	    Fy = rho*X'*Py*X
	    Fz = rho*X'*Pz*X

	    Mx = rho*X'*Qx*X
	    My = rho*X'*Qy*X
	    Mz = rho*X'*Qz*X


	    # state derivative calculations
	    #dxst = zeros(12)
	    dxst[1] = r*v-q*w+Fx/m-g*sin(theta)+T/m
	    dxst[2] = p*w-r*u+Fy/m+g*cos(theta)*sin(phi)  # Force Equations
	    dxst[3] = q*u-p*v+Fz/m+g*cos(theta)*cos(phi)

	    dxst[4] = (c1*r+c2*p-c4*Ip*OmegaP)*q+c3*Mx+c4*Mz
	    dxst[5] = (c5*p+c7*Ip*OmegaP)*r-c6*(p^2-r^2)+c7*My   # moment equations
	    dxst[6] = (c8*p-c2*r-c9*Ip*OmegaP)*q+c9*Mz+c4*Mx

	    dxst[7] = p+tan(theta)*(q*sin(phi)+r*cos(phi))
	    dxst[8] = q*cos(phi)-r*sin(phi)
	    dxst[9] = (q*sin(phi)+r*cos(phi))/cos(theta)

	    dxst[10:12] = [1 0 0; 0 cos(phi) sin(phi);0 -sin(phi) cos(phi)]*[cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)]*[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0; 0 0 1]*X[1:3]


	end


	prob = ODEProblem(calc_dx!,xst0,tspan,pq)

	#nsteps = 111
	#tstops = LinRange(tspan[1],tspan[2],nsteps)
	#sol = solve(prob,tstops=tstops,adaptive=false)
	#sol  = solve(prob,Rodas4())
	#sol  = solve(prob,Rosenbrock23())
	sol = solve(prob,dtmin=1e-6,Rosenbrock23(),saveat=LinRange(t_begin,t_end,tsteps))
	return sol
end
"""
    main(inp_file)

	Main function to calculate all the values required

	Input : inp_file - Input File

	Output:
"""
function main(inp_file,iseq=0)

    sref,bref,cref,aseq,vinf,ρ,alpha,nwing,
	ntot_patch,npatch_perwing,
	ntot_lat,nspan_perpatch,nlat_perpatch,xyz_qc_patch,
	chord_patch,twist_patch,α_zl_patch = read_inp(inp_file)

	nspn,spn_map,spn_loc,sbar,ebar,mbar,nbar,cbar,dbar,tbar,chord,twist,α_zl,θ,Λ,ds = geom_calc(nwing,ntot_patch,npatch_perwing,ntot_lat,nspan_perpatch,nlat_perpatch,xyz_qc_patch,chord_patch,twist_patch,α_zl_patch)

	AIC,AICₘ,AICₜ,J,Px,Py,Pz,Qx,Qy,Qz,Gxgamma,Gygamma,
	Gzgamma,Hxgamma,Hygamma,Hzgamma,Gxv,Gyv,Gzv,Hxv,Hyv,Hzv,Gxomega,
	Gyomega,Gzomega,Hxomega,Hyomega,Hzomega,crsbar = calc_AICs(ntot_lat,sbar,ebar,cbar,mbar,nbar,tbar)

	sol = calc_motion(Px,Py,Pz,Qx,Qy,Qz)

	return sol,Px,Py,Pz,Qx,Qy,Qz
end



#Cl_spn,nspn,spn_map,spn_loc,θ,rhs,AIC,AICₜ,AICₘ,Λ,sbar,ebar,ds,Γ,chord,cbar,mbar,nbar,Cl,CL,CDind,CDind_ff,SpnLd = main("input.yaml");
#Px = Py = Pz = Qx = Qy = Qz = zeros(6,6)
sol,Px,Py,Pz,Qx,Qy,Qz = main("input.yaml");

nsteps= size(sol.u,1)
xyz = zeros(3,nsteps)

for i in 1:nsteps
    xyz[:,i] = sol.u[i][10:12]
end


plt = plot3d(1,xlabel="X-coord",ylabel="Y-coord",zlabel="Z-coord",xlim=(minimum(xyz[1,:]),maximum(xyz[1,:])),ylim=(minimum(xyz[2,:]),maximum(xyz[2,:])),zlim=(minimum(xyz[3,:]),maximum(xyz[3,:])),marker=2)
@gif for i in 1:nsteps
    push!(plt,xyz[1,i],xyz[2,i],xyz[3,i])
end every 1
