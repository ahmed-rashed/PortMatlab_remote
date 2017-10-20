MODULE MEX_INTERFACES
INTERFACE
	subroutine mexErrMsgTxt(errormsg)
		character*(*) errormsg
	end subroutine mexErrMsgTxt

	subroutine mexWarnMsgTxt(warningmsg)
		character*(*) warningmsg
	end subroutine mexWarnMsgTxt


END INTERFACE
END MODULE MEX_INTERFACES