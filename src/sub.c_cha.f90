! Y. Yagi
character function cha(nm)
  integer,intent(in)  :: nm
  if(nm.eq.0) cha='0'
  if(nm.eq.1) cha='1'
  if(nm.eq.2) cha='2'
  if(nm.eq.3) cha='3'
  if(nm.eq.4) cha='4'
  if(nm.eq.5) cha='5'
  if(nm.eq.6) cha='6'
  if(nm.eq.7) cha='7'
  if(nm.eq.8) cha='8'
  if(nm.eq.9) cha='9'
  return
END function cha

