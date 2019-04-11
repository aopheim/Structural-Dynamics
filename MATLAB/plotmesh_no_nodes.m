function  outputstring=plotmesh_no_nodes(Nodes,Elements),

Nel = size(Elements,1);

clf;        %Clears the current figure
%plot3(Nodes(:,1),Nodes(:,2),Nodes(:,3),'k*')
for el=1:Nel,
   type   = Elements(el,1);
   node1  = Nodes(Elements(el,3),:);
   node2  = Nodes(Elements(el,4),:);
   hold on
   if type==1,              % if bar
      plot3([node1(1);node2(1)],[node1(2);node2(2)],[node1(3);node2(3)],':k', 'Linewidth', 2);
   elseif type==2,          %if beam
      plot3([node1(1);node2(1)],[node1(2);node2(2)],[node1(3);node2(3)],':k', 'Linewidth', 2);
   end
end
axis equal
%axis off
view(-65,35)
