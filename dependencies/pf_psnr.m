function[ratio] = pf_psnr(y, y_tilde)
	m = length(y);

	if(mod(m, 3) == 0)
		seq1 = 1:(m/3); seq2 = ((m/3)+1):(2*(m/3)); seq3 = (2*(m/3)+1):m;
		max1 = max(y_tilde(seq1)); max2 = max(y_tilde(seq2)); max3 = max(y_tilde(seq3));
		
		psnr1 = psnr(y(seq1), y_tilde(seq1), max1);
		psnr2 = psnr(y(seq2), y_tilde(seq2), max2);
		psnr3 = psnr(y(seq3), y_tilde(seq3), max3);
		
		ratio = (psnr1+psnr2+psnr3)/3;

	elseif(mod(m, 4) == 0)
		seq1 = 1:(m/4); seq2 = ((m/4)+1):(2*(m/4)); 
		seq3 = (2*(m/4)+1):(3*(m/4)); seq4 = (3*(m/4)+1):m;
		
		max1 = max(y_tilde(seq1)); max2 = max(y_tilde(seq2)); 
		max3 = max(y_tilde(seq3)); max4 = max(y_tilde(seq4));
		
		psnr1 = psnr(y(seq1), y_tilde(seq1), max1);
		psnr2 = psnr(y(seq2), y_tilde(seq2), max2);
		psnr3 = psnr(y(seq3), y_tilde(seq3), max3);
		psnr4 = psnr(y(seq4), y_tilde(seq4), max4);

		ratio = (psnr1+psnr2+psnr3+psnr4)/4;
		
	else 
		ratio = psnr(y, y_tilde, max(y_tilde));
	end

end