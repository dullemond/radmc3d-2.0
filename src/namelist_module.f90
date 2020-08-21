module namelist_module

contains

!===================================================================
!                   ROUTINES FOR PARAMETER LIST READING
!===================================================================


!-------------------------------------------------------------------
!            READ IN THE ENTIRE INPUT FILE INTO A STRING ARRAY
!-------------------------------------------------------------------
subroutine read_input_file()
  implicit none
  !
  common/parses/inpstring
  save /parses/
  character*160 inpstring(300)
  common/parsei/nrstring,lineok
  save /parsei/
  
  integer nrstring,lineok(300)
  !
  character*160 string
  integer iline
  !
  do iline=1,300
     read(1,'(A158)',end=10) string
     inpstring(iline) = string
     lineok(iline) = 0
  enddo
  write(*,*) 'Input file contains too many lines...'
  stop
10 continue
  nrstring=iline-1
  !      
end subroutine read_input_file



!-------------------------------------------------------------------
!         BASIC PARSING ROUTINE: SPLIT STRING INTO NAME AND VALUE
!
!     This routine interprets lines such as
!  
!       ; The following variable is useful:
!       var1   = 5.74d0  ; Here some comments
!
!     The first line is simply ignored (though marked as OK) and
!     the second line would be split into strname='var1' and
!     strvalue='5.74d0' with lenname being 4 and lenvalue being 6!
!     
!-------------------------------------------------------------------
subroutine parse_name_value(string,strname,lenname,strvalue,lenvalue,icomm)
  implicit none
  character*160 string
  character*80 strname,strvalue
  integer lenname,lenvalue
  logical icomm
  !
  integer ilen,i,istart,istop
  !
  icomm    = .false.
  strname  = ''
  lenname  = 0
  strvalue = ''
  lenvalue = 0
  ilen     = len_trim(string)
  istart   = 1
  do i=1,ilen-1
     if(string(i:i) == ' ') then
        istart=i+1
     else
        goto 20
     endif
  enddo
  icomm=.true.
  return
20 continue
  if(string(istart:istart) == ';') then
     icomm=.true.
     return
  endif
  if(istart >= ilen) return
  do i=istart+1,ilen-1
     if(string(i:i) == ' ') goto 40
     if(string(i:i) == '=') goto 40
  enddo
  return
40 continue
  istop=i-1
  lenname = istop-istart+1
  if(lenname > 80) then
     write(*,*) 'Variable name too long'
     stop
  endif
  strname = string(istart:istop)
  do i=istop+1,ilen-1
     if(string(i:i) == '=') goto 50
  enddo
  lenname=0
  strname=''
  lenvalue=0
  strvalue='' 
  return
50 continue
  istart = i+1
  if(istart > ilen) stop
  do i=istart,ilen
     if(string(i:i).eq.' ') then
        istart=i+1
     else
        goto 60
     endif
  enddo
60 continue
  do i=istart+1,ilen
     if(string(i:i) == ';') goto 80
     if(string(i:i) == ' ') goto 80
  enddo
80 continue
  istop = i-1
  lenvalue = istop-istart+1
  if(lenvalue == 80) then
     write(*,*) 'Value length too long'
     stop
  endif
  strvalue = string(istart:istop)
  !
  !      write(*,*) '##',strname(1:lenname),'##',strvalue(1:lenvalue),'##'
  !     
end subroutine parse_name_value


!-------------------------------------------------------------------
!                         STRING COMPARING
!-------------------------------------------------------------------
function stringcompare(string1,string2,len)
  implicit none
  character*80 :: string1,string2
  integer :: len
  integer :: i
  logical :: stringcompare
  !
  if(len > 80) stop
  if(len <= 0) then
     stringcompare=.false.
     return
  endif
  stringcompare=.true.
  do i=1,len
     if(string1(i:i) /= string2(i:i)) stringcompare=.false.
  enddo
  return
end function stringcompare


!-------------------------------------------------------------------
!                       PARSING ROUTINE FOR DOUBLE
!-------------------------------------------------------------------
subroutine parse_input_double(nameorig,value)
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
  double precision :: value
  character*30 nameorig
  character*80 name
  !
  common/parses/inpstring
  save /parses/
  character*160 inpstring(300)
  common/parsei/nrstring,lineok
  save /parsei/
  integer nrstring,lineok(300)
  !
  integer iline,ichar
  character*80 strname,strvalue
  integer lenname,lenvalue,lenn
  logical found,icomm
  !
  name(1:30)=nameorig
  found=.false.
  do ichar=2,30
     if(name(ichar:ichar) == '@') goto 30
  enddo
  write(*,*) 'INTERNAL ERROR: Keywords must end with @!!!'
  stop
30 continue
  lenn=ichar-1
  do iline=1,nrstring
     call parse_name_value(inpstring(iline), strname,lenname,           &
          strvalue, lenvalue,icomm)
     if(icomm) then
        lineok(iline) = 1
     elseif(stringcompare(strname,name,lenn) .and.                      &
          (lenn.eq.lenname)) then
        if(found) then
           write(*,*) 'Found one parameter more than once'
           write(*,*) 'in the input file: ',strname(1:lenname)
           write(*,*) '*** ABORTING ***'
           stop
        endif
        read(strvalue(1:lenvalue),*) value
        found=.true.
        lineok(iline) = 1
     endif
  enddo
  !
end subroutine parse_input_double

!-------------------------------------------------------------------
!                       PARSING ROUTINE FOR INTEGER
!-------------------------------------------------------------------
subroutine parse_input_integer(nameorig,value)
  implicit none
  character*30 nameorig
  character*80 name
  integer value
  !
  common/parses/inpstring
  save/parses/
  character*160 inpstring(300)
  common/parsei/nrstring,lineok
  save/parsei/
  integer nrstring,lineok(300)
  !
  integer iline,ichar
  character*80 strname,strvalue
  integer lenname,lenvalue,lenn
  logical found,icomm
  !
  name(1:30)=nameorig
  found=.false.
  do ichar=2,30
     if(name(ichar:ichar) == '@') goto 30
  enddo
  write(*,*) 'INTERNAL ERROR: Keywords must end with @!!!'
  stop
30 continue
  lenn=ichar-1
  do iline=1,nrstring
     call parse_name_value(inpstring(iline), strname,lenname,           &
          strvalue,lenvalue,icomm)
     if(icomm) then
        lineok(iline) = 1
     elseif(stringcompare(strname,name,lenn) .and.                      &
          (lenn == lenname)) then
        if(found) then
           write(*,*) 'Found one parameter more than once'
           write(*,*) 'in the input file: ',strname(1:lenname)
           write(*,*) '*** ABORTING ***'
           stop
        endif
        read(strvalue(1:lenvalue),*) value
        found=.true.
        lineok(iline) = 1
     endif
  enddo
  !
end subroutine parse_input_integer

!-------------------------------------------------------------------
!                       PARSING ROUTINE FOR INTEGER*8
!-------------------------------------------------------------------
subroutine parse_input_integer8(nameorig,value)
  implicit none
  character*30 nameorig
  character*80 name
  integer*8 :: value
  !
  common/parses/inpstring
  save/parses/
  character*160 inpstring(300)
  common/parsei/nrstring,lineok
  save/parsei/
  integer nrstring,lineok(300)
  !
  integer iline,ichar
  character*80 strname,strvalue
  integer lenname,lenvalue,lenn
  logical found,icomm
  !
  name(1:30)=nameorig
  found=.false.
  do ichar=2,30
     if(name(ichar:ichar) == '@') goto 30
  enddo
  write(*,*) 'INTERNAL ERROR: Keywords must end with @!!!'
  stop
30 continue
  lenn=ichar-1
  do iline=1,nrstring
     call parse_name_value(inpstring(iline), strname,lenname,           &
          strvalue,lenvalue,icomm)
     if(icomm) then
        lineok(iline) = 1
     elseif(stringcompare(strname,name,lenn) .and.                      &
          (lenn == lenname)) then
        if(found) then
           write(*,*) 'Found one parameter more than once'
           write(*,*) 'in the input file: ',strname(1:lenname)
           write(*,*) '*** ABORTING ***'
           stop
        endif
        read(strvalue(1:lenvalue),*) value
        found=.true.
        lineok(iline) = 1
     endif
  enddo
  !
end subroutine parse_input_integer8

!-------------------------------------------------------------------
!                       PARSING ROUTINE FOR WORD
!-------------------------------------------------------------------
subroutine parse_input_word(nameorig,value,len)
  implicit none
  character*30 nameorig
  character*80 name
  character*90 value
  integer len
  !
  common/parses/inpstring
  save/parses/
  character*160 inpstring(300)
  common/parsei/nrstring,lineok
  save/parsei/
  integer nrstring,lineok(300)
  !
  integer iline,ichar
  character*80 strname,strvalue
  integer lenname,lenvalue,lenn
  logical found,icomm
  !
  name(1:30)=nameorig
  found=.false.
  do ichar=2,30
     if(name(ichar:ichar) == '@') goto 30
  enddo
  write(*,*) 'INTERNAL ERROR: Keywords must end with @!!!'
  stop
30 continue
  lenn=ichar-1
  do iline=1,nrstring
     call parse_name_value(inpstring(iline), strname,lenname,           &
          strvalue,lenvalue,icomm)
     if(icomm) then
        lineok(iline) = 1
     elseif(stringcompare(strname,name,lenn) .and.                      &
          (lenn == lenname)) then
        if(found) then
           write(*,*) 'Found one parameter more than once'
           write(*,*) 'in the input file: ',strname(1:lenname)
           write(*,*) '*** ABORTING ***'
           stop
        endif
        read(strvalue(1:lenvalue),*) value
        len = lenvalue
        found=.true.
        lineok(iline) = 1
     endif
  enddo
  !
end subroutine parse_input_word

!-------------------------------------------------------------------
!                   CHECK IF ALL LINES ARE OK
!-------------------------------------------------------------------
subroutine check_all_lines_ok(ierror)
  implicit none
  logical ierror
  !
  common/parses/inpstring
  save/parses/
  character*160 inpstring(300)
  common/parsei/nrstring,lineok
  save/parsei/
  integer nrstring,lineok(300)
  !
  character*160 string
  integer iline,ilen
  !
  ierror = .false.
  do iline=1,nrstring
     if(lineok(iline).eq.0) then
        write(*,*) 'ERROR in input file: variable unknown:'
        ilen = len_trim(inpstring(iline))
        string = inpstring(iline)
        write(*,*) string(1:ilen)
        ierror = .true.
     endif
  enddo
  !
end subroutine check_all_lines_ok


end module namelist_module
