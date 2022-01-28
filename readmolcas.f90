program main
  implicit none
  character(len=25) :: radiusfile, energyfile, dipolefile, polfile, rotconstfile,polxyzfile
  character(len=100) :: outputfile, pesfile
  character(len=100) :: comando, auxchprev
  character(len=:), allocatable :: auxch
  integer :: lineas, palabras, caracteres
  integer :: ij, jk, ik, fin
  logical :: existe
  integer :: lineaspol
  real(kind(1.0d0)), allocatable :: r(:), energy(:)
  real(kind(1.0d0)), allocatable :: dipole(:,:) !! #1: 1 corresponds to x, 2 to y, 3 to z and 4 to total.
  real(kind(1.0d0)), allocatable :: polarizabilities(:,:) !! diagonal elements of the polarizability
  !! #1: 1 corresponds to xx, 2 to yy, 3 to zz and 4 to total.
  real(kind(1.0d0)), parameter :: zero=0
  real(kind(1.0d0)), allocatable :: alpha(:,:) !! #1 and #2: 1 corresponds to x, 2 to y and 3 to z
  real(kind(1.0d0)), allocatable :: rotconst(:,:) !! Rotational constants


  !! Initialize
  lineas=0
  palabras=0
  caracteres=0
  write(outputfile,'(A)') 'output.out'
  write(radiusfile,'(A)') 'radius.txt'
  write(energyfile,'(A)') 'energy.txt'
  write(dipolefile,'(A)') 'dipole.txt'
  write(polfile,'(A)') 'polarizabilities.txt'
  write(polxyzfile,'(A)') 'alphaxyz.dat'
  write(rotconstfile,'(A)') 'rotconstfile.txt'

  !! Remove auxiliar files
  
  write(comando,'(A,A)') "rm -f ", radiusfile

  CALL EXECUTE_COMMAND_LINE(comando)
  
  write(comando,'(A,A)') "rm -f ", energyfile

  CALL EXECUTE_COMMAND_LINE(comando)
  
  write(comando,'(A,A)') "rm -f ", dipolefile

  CALL EXECUTE_COMMAND_LINE(comando)

  !! Call the input file

  call getarg(1,outputfile)
  INQUIRE(FILE=trim(outputfile), EXIST=existe)
  if (.not.existe) then
     write(*,*) 'The file ', trim(outputfile), ' does not exist.'
     stop
  end if
  
  !! Store the radius, energy and dipole files
  
  write(comando,'(4A)') "grep -i 'EXPORT' ", trim(outputfile), " >> ", radiusfile

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

  write(comando, '(A,A,A)') "wc ", radiusfile, " > auxiliarfile" 
  
  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

  open(123, file='auxiliarfile')

  read(123, *) lineas, palabras, caracteres, auxchprev

  close(123)

  !! now read the energies

  write(comando,'(A,A)') "rm -f ", energyfile

  CALL EXECUTE_COMMAND_LINE(comando)
  
  write(comando,'(4A)') "grep -i 'energy=' ", trim(outputfile), " >> ", energyfile

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

  !! now read the dipoles

  write(comando,'(A,A)') "rm -f ", dipolefile

  CALL EXECUTE_COMMAND_LINE(comando)
  
  write(comando,'(4A)') "grep -i 'Total=' ", trim(outputfile), " >> ", dipolefile

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

  !! now read the polarizabilities

  write(comando,'(A,A)') "rm -f ", polfile

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

  write(comando,'(4A)') "grep -i -n 'Polarizabilities' ", trim(outputfile), " >> ", polfile

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

  write(comando,'(3A)') "wc ", trim(polfile), " >> fort.111"

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

  !! if the polfile is empty then we delete it
  
  ij=1
  
  read(111,*) ij

  if (ij.eq.0) then
  write(comando,'(2A)') "rm -rf ", polfile

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

  write(comando,'(2A)') "rm -rf ", rotconstfile

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)
     
  end if
  
  !! now read the rotational constants

  write(comando,'(A,A)') "rm -f ", rotconstfile

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)
  
  write(comando,'(4A)') "grep -i 'Rotational Constants (GHz)' ", trim(outputfile), " >> ", rotconstfile

  CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

  !! The number of distances computed is "lineas". Now we set an array to store the distances, energies
  !! dipole moments and polarizabilities.

  !! Store
  allocate(r(1:lineas))
  allocate(energy(1:lineas))
  allocate(dipole(1:4,1:lineas))
  allocate(polarizabilities(1:3,1:lineas))
  allocate(alpha(1:3, 1:3))
  allocate(rotconst(1:3,1:lineas))
  
  !! initialize
  r=zero
  energy=zero
  dipole=zero
  polarizabilities=zero
  rotconst=zero
  
  open(124,file=radiusfile,status='old')
  open(125,file=energyfile,status='old')
  open(126,file=dipolefile,status='old')

  !! check if the polarizability file exist
  
  Inquire(FILE=trim(polfile), EXIST=existe)
  if (existe) open(127,file=polfile,status='old')
  
  open(129,file=rotconstfile,status='old')

  do jk=1,lineas

     write(*,'(A,I2.1, A, I2.2)') 'Reading coordinate ', jk, ' out of ', lineas
     
     read(124,'(A)') auxchprev

     auxch=trim(auxchprev)

     !! store the radius
     Do ij=len(auxch),1,-1
        !! search for the number 
        if (auxch(ij:ij).eq.'.') then
           !! found, store it in the auxiliar variable
           auxchprev=trim(auxch(ij-1:len(auxch)))
           !! save it in the array
           read(auxchprev,*) r(jk)
           exit
        end if
     End Do

     !! store the energy
     read(125,*) auxchprev,energy(jk)
          
     !! store the dipole
     !! the first one corresponds to SCF
     read(126,*) auxchprev,dipole(1,jk),auxchprev,dipole(2,jk),auxchprev,dipole(3,jk),auxchprev,dipole(4,jk)
     !! this one corresponds to RASSCF if we only compute one orbital
     read(126,*) auxchprev,dipole(1,jk),auxchprev,dipole(2,jk),auxchprev,dipole(3,jk),auxchprev,dipole(4,jk)

     !! look for the polarizabilities
     !! first we look for alphaxx

     !! we read line where the polarizability section starts

     Inquire(FILE=trim(polfile), EXIST=existe)
     if (existe) then
     read(127,*,IOSTAT=fin) auxchprev

     if (fin.ne.0) then
        write(*,*) 'End of the polarizability file reached'
        exit !! Exit the loop in jk to runs over the file
     end if
     
     auxch=trim(auxchprev)
     Do ij=1,len(auxch)
        !! search for the number 
        if (auxch(ij:ij).eq.':') then
           !! found, store it in the auxiliar variable
           auxchprev=trim(auxch(1:ij-1))
           !! save it in the array
           read(auxchprev,*) lineaspol
           exit
        end if
     End Do


     !!write the polarizability tensor on a file
     !!first row
     write(comando,*) "head -n ", lineaspol+6, trim(outputfile), "| tail -1 >", polxyzfile
     CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

     !!second row
     write(comando,*) "head -n ", lineaspol+7, trim(outputfile), "| tail -1 >>", polxyzfile
     CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)

     !!third row
     write(comando,*) "head -n ", lineaspol+8, trim(outputfile), "| tail -1 >>", polxyzfile
     CALL EXECUTE_COMMAND_LINE(comando,wait=.true.)
       
     alpha=zero

     !!store the polarizability
     open(128,file=polxyzfile,status='old')

     read(128,*) alpha(1,1)
     read(128,*) alpha(2,1), alpha(2,2)
     read(128,*) alpha(3,1), alpha(3,2), alpha(3,3)
     
     close(128)

     do ik=1,3
        polarizabilities(ik,jk)=alpha(ik,ik)
     end do     
     
     !!now we look for the rotational constants and store them in the array

     rotconst(:,jk)=zero
     read(129,'(A31,2E8.4)') auxchprev, (rotconst(ij,jk), ij=1,2)
  end if

end do
  
  !! Now we write everything in a file
!  close(124)
!  close(125)
!  close(126)
!  close(127)
!  close(129)

  !! set output file

  write(pesfile,'(A4,A)') 'PES_', trim(outputfile)

  Open(130, file=trim(pesfile))

  write(130,'(12A16)') '# R', 'E (hartree)', 'mu_x (D)','mu_y (D)','mu_z (D)','mu_t (D)', 'alphaxx','alphayy','alphazz', &
       'Bxx (GHz)', 'Byy (GHz)'
  
  Do jk=1,lineas
     write(130,'(12ES16.6E2)') r(jk), energy(jk), (dipole(ij,jk), ij=1,4), &
          (polarizabilities(ij,jk), ij=1,3), (rotconst(ij,jk), ij=1,2)
  End Do

  close(130)
  
end program main
