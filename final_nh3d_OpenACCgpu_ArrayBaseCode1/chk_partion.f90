!!  A test code for processor partions ith 2d MPI codes
 Program Chk_Partition
 Implicit None 
 
   Integer, Parameter :: nproc=8 , nex=36, ney=16 
   Integer, Parameter :: blocksx=4 , blocksy=2
   integer :: part(nproc) 
   integer :: rank, r,k, ix,iy, colx,coly,extra,offsetx  ,offsety
   Integer ::    rankS, rankN , rankW, rankE



    do r = 0, nproc-1 

      ix = modulo(r,blocksx)
      iy = rank / blocksx ! integer div intended

        colx = nex / blocksx ! integer div intentional
        extra = modulo(nex, blocksx)
        if (ix < extra) then
          colx = colx + 1
          offsetx = colx * ix
        else
          offsetx = extra + colx * ix
        end if

        ! blocks in y
        coly = ney / blocksy
        extra = modulo(ney, blocksy)
        if (iy < extra) then
          coly = coly + 1
          offsety = coly * iy
        else
          offsety = extra + coly * iy
        end if

        rank =r 

        rankS = modulo(rank - blocksx, nproc)
        rankE = iy*blocksx + modulo(rank + 1, blocksx)
        rankN = modulo(rank + blocksx, nproc)
        rankW = iy*blocksx + modulo(rank - 1, blocksx)

      ! print*, rank, ix, iy 
        print('(a,i3,a,i3,i3,a,i5,i5)'),'rank :', rank,", colmn x/y:", colx, coly , &
                        ", offset x/y:" ,offsetx , offsety 

      ! Processor locality for torus topology
      ! print*, 'S&N',  rankS, rankN
      ! print*, 'W&E',  rankW, rankE
    end do


      
      !  do r = 0, 7
      !   k = modulo(r-4,8)
      !   print*, r, k , floor(dble(r)/4.0)
      !  enddo 
 end program Chk_Partition  
