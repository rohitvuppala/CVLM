include("functions_cvlm.jl")

#Cl_spn,nspn,spn_map,spn_loc,θ,rhs,AIC,AICₜ,AICₘ,Λ,sbar,ebar,ds,Γ,chord,cbar,mbar,nbar,Cl,CL,CDind,CDind_ff,SpnLd = main("input.yaml");
#Px = Py = Pz = Qx = Qy = Qz = zeros(6,6)
sol,Px,Py,Pz,Qx,Qy,Qz,lat2wing_map = main("input.yaml","motion_param.yaml");

nsteps= size(sol.u,1)
xyz = zeros(3,nsteps)
xyz1= zeros(3,nsteps)
for i in 1:nsteps
    xyz[:,i]  = sol.u[i][10:12]
    xyz1[:,i] = solwing.u[i][10:12]
end


p = plot3d(2,labels=["Wing+tail" "Wing only"],xlabel="X-coord",ylabel="Y-coord",zlabel="Z-coord",xlim=(minimum(xyz[1,:]),maximum(xyz[1,:])),ylim=(minimum(xyz[2,:]),maximum(xyz[2,:])),zlim=(minimum(xyz[3,:]),maximum(xyz[3,:])),marker=2)
#p = plot3d!(p,2)
@gif for i in 1:nsteps
    push!(p,1,xyz[1,i],xyz[2,i],xyz[3,i])
    push!(p,2,xyz1[1,i],xyz1[2,i],xyz1[3,i])
end every 1