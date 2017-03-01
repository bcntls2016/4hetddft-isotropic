 function Morse_HeTiO2_1D(r) result( output )
	implicit none
! 	real(8), intent(in) :: x
! 	real(8), intent(in) :: y
! 	real(8), intent(in) :: z
	real(8) :: output
	
	real(8), parameter :: De = 39.95_8
	real(8), parameter :: alpha = 1.721_8
	real(8), parameter :: Re = 4.3619_8
	
	real(8), intent(in) :: r
	
! 	r = sqrt(sum([x,y,z]**2))
	output = De*( exp(2.0_8*alpha*(Re-r)) - 2.0_8*exp(alpha*(Re-r)) )*1.4387752
end function Morse_HeTiO2_1D


!  function Morse_HeTiO2( x, y, z ) result( output )
! 	implicit none
! 	real(8), intent(in) :: x
! 	real(8), intent(in) :: y
! 	real(8), intent(in) :: z
! 	real(8) :: output
! 	
! 	real(8), parameter :: De = 39.95_8
! 	real(8), parameter :: alpha = 1.721_8
! 	real(8), parameter :: Re = 4.3619_8
! 	
! 	real(8) :: r
! 	
! 	r = sqrt(sum([x,y,z]**2))
! 	output = De*( exp(2.0_8*alpha*(Re-r)) - 2.0_8*exp(alpha*(Re-r)) )
! end function Morse_HeTiO2

	function Morse_HeTiO2_3D( x, y, z ) result( output )
		implicit none
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output

		integer :: i,j

		real(8) :: pi
		real(8) :: a, b
		real(8) :: De(3,4)
		real(8) :: alpha(3,4)
		real(8) :: Ze(3,4)

		real(8) :: Lx(3,1)
		real(8) :: Ly(4,1)
		real(8) :: De_xy(1,1)
		real(8) :: alpha_xy(1,1)
		real(8) :: Ze_xy(1,1)

		real(8) :: rbuffer(1,1)

		pi = acos(-1.0_8)

		a = 2.9850_8
		b = 6.49548_8

		De = transpose(reshape( &
			[ 54.671500_8, -2.228480_8, 0.585194_8, 13.211000_8, &
			   4.842040_8,  0.363525_8, 4.614530_8, -1.182030_8, &
			   3.747220_8,  3.222690_8, 0.333350_8,  0.587320_8  ], &
		     [4,3] ))

		alpha = transpose(reshape( &
			[  1.492420_8, -0.010921_8, 0.003189_8, -0.168927_8, &
			  -0.097618_8, -0.024780_8, 0.062441_8,  0.016801_8, &
			  -0.034535_8,  0.030015_8, 0.012846_8, -0.011699_8  ], &
		     [4,3] ))

		Ze = transpose(reshape( &
			[  4.073990_8,  0.057537_8, -0.001500_8, -0.568990_8, &
			  -0.115578_8, -0.018127_8, -0.090046_8,  0.008590_8, &
			  -0.075153_8, -0.051210_8, -0.023349_8, -0.017931_8  ], &
		     [4,3] ))
		
		Lx = transpose(reshape( &
			cos( 2*pi*[ 0.0_8, 1.0_8, 2.0_8 ]*x/a ), &
		     [1,3] ))

		Ly = transpose(reshape( &
			cos( 2*pi*[ 0.0_8, 1.0_8, 2.0_8, 3.0_8 ]*y/b ), &
		     [1,4] ))

		De_xy = matmul( matmul( transpose(Lx), De ), Ly )
		alpha_xy = matmul( matmul( transpose(Lx), alpha ), Ly )
		Ze_xy = matmul( matmul( transpose(Lx), Ze ), Ly )
	
		rbuffer = De_xy*( exp(-2.0_8*alpha_xy*( z-Ze_xy ))-2.0_8*exp(-alpha_xy*(z-Ze_xy)) )

		output = rbuffer(1,1)
	end function Morse_HeTiO2_3D